/** Orekit Code for Orbit Propagator
 * Add description 
 * 
 */
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.hipparchus.ode.AbstractIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.ForceModel;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.conversion.FiniteDifferencePropagatorConverter;
import org.orekit.propagation.conversion.KeplerianPropagatorBuilder;
import org.orekit.propagation.conversion.PropagatorBuilder;
import org.orekit.propagation.conversion.PropagatorConverter;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

public class OrbitPropagator {
	
    public static void main(final String[] args) {
        try {
        	
            // configure Orekit
        	File orekitData = new File("C:\\Users\\Smit\\Desktop\\Orekit\\orekit-data");
        	DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
        	manager.addProvider(new DirectoryCrawler(orekitData));
            final File home       = new File(System.getProperty("user.home"));
            if (!orekitData.exists()) {
                System.err.format(Locale.US, "Failed to find %s folder%n",
                                  orekitData.getAbsolutePath());
                System.err.format(Locale.US, "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n",
                                  "orekit-data-master.zip", "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip",
                                  home.getAbsolutePath());
                System.exit(1);
            }
            manager.addProvider(new DirectoryCrawler(orekitData));

            // gravity field
            final NormalizedSphericalHarmonicsProvider provider = GravityFieldFactory.getNormalizedProvider(2, 0);
            final double mu =  provider.getMu();

            // inertial frame
            final Frame inertialFrame = FramesFactory.getGCRF();

            // Initial date
            final AbsoluteDate initialDate = new AbsoluteDate(2019, 8, 22, 0, 0, 00.000,
                                                              TimeScalesFactory.getUTC());

            // Initial orbit (GTO)
            final double a     = 7100000;                // semi major axis in meters
            final double e     = 0.05;              // eccentricity
            final double i     = FastMath.toRadians(30);   // inclination
            final double omega = FastMath.toRadians(0.0001); // perigee argument
            final double raan  = FastMath.toRadians(50); // right ascention of ascending node
            final double lM    = 0;                       // true anomaly
            final double muc = 3.986004418e14;
            final Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.TRUE,
                                                          inertialFrame, initialDate, mu);
            final double period = initialOrbit.getKeplerianPeriod();
            
            // Initial state definition
            final SpacecraftState initialState = new SpacecraftState(initialOrbit);

            // Adaptive step integrator with a minimum step of 0.001 and a maximum step of 1000
            final double minStep    = 0.001;
            final double maxStep    = 1000.;
            final double dP         = 1.e-2;
            final double dura = 172800 ;
            final OrbitType orbType = OrbitType.CARTESIAN;
            final double[][] tol = NumericalPropagator.tolerances(dP, initialOrbit, orbType);
            final AbstractIntegrator integrator = new DormandPrince853Integrator(minStep, maxStep,
                                                                                 tol[0], tol[1]);

            // Propagator
            final NumericalPropagator numProp = new NumericalPropagator(integrator);
            numProp.setInitialState(initialState);
            numProp.setOrbitType(orbType);

            // Force Models:
            // ( only J2 is considered here)
            // final ForceModel gravity = new HolmesFeatherstoneAttractionModel(FramesFactory.getITRF(IERSConventions.IERS_2010, true), provider);
            final ForceModel gravity = new NewtonianAttraction( muc);
            
            // Add force models to the propagator
            numProp.addForceModel(gravity);

            // Propagator factory
            final PropagatorBuilder builder = new KeplerianPropagatorBuilder(initialOrbit, PositionAngle.TRUE, dP);

            // Propagator converter
            final PropagatorConverter fitter = new FiniteDifferencePropagatorConverter(builder, 1.e-6, 5000);

            // Resulting propagator
            final KeplerianPropagator kepProp = (KeplerianPropagator) fitter.convert(numProp, 2 *  period, 251);

            // Step handlers
            final StatesHandler numStepHandler = new StatesHandler();
            final StatesHandler kepStepHandler = new StatesHandler();

            // Set up fixed step handlers to the propagator
            numProp.getMultiplexer().add(60., numStepHandler);
            kepProp.getMultiplexer().add(60., kepStepHandler);

            // Extrapolate from the initial to the final date
            numProp.propagate(initialDate, initialDate.shiftedBy(dura));
            kepProp.propagate(initialDate, initialDate.shiftedBy(dura));

            // retrieve the states
            final List<SpacecraftState> numStates = numStepHandler.getStates();
            final List<SpacecraftState> kepStates = kepStepHandler.getStates();

            // Print the results on the output file
            final File output = new File(home, "elements.txt");
            try (PrintStream stream = new PrintStream(output, StandardCharsets.UTF_8.name())) {
                // stream.println("# date Anum Akep Enum Ekep Inum Ikep LMnum LMkep");
                for (SpacecraftState numState : numStates) {
                    for (SpacecraftState kepState : kepStates) {
                        if (numState.getDate().compareTo(kepState.getDate()) == 0) { 
                            stream.println( // numState.getDate() +
                                           " " + numState.getA() +
                                           " " + kepState.getA() +
                                           " " + numState.getE() +
                                           " " + kepState.getE() +
                                           " " + FastMath.toDegrees(numState.getI()) +
                                           " " + FastMath.toDegrees(kepState.getI()) +
                                           " " + FastMath.toDegrees(MathUtils.normalizeAngle(numState.getLM(), FastMath.PI)) +
                                           " " + FastMath.toDegrees(MathUtils.normalizeAngle(kepState.getLM(), FastMath.PI)));
                            break;
                        }
                    }
                }
            }
            final String resultSavedAs = "Results saved as file ";
            System.out.println(resultSavedAs + output);

            final File output1 = new File(home, "elts_pv.txt");
            try (PrintStream stream = new PrintStream(output1, StandardCharsets.UTF_8.name())) {
                // stream.println("# date pxn pyn pzn vxn vyn vzn pxk pyk pzk vxk vyk vzk");
                for (SpacecraftState numState : numStates) {
                    for (SpacecraftState kepState : kepStates) {
                        if (numState.getDate().compareTo(kepState.getDate()) == 0) {
                            final double pxn = numState.getPVCoordinates().getPosition().getX();
                            final double pyn = numState.getPVCoordinates().getPosition().getY();
                            final double pzn = numState.getPVCoordinates().getPosition().getZ();
                            final double vxn = numState.getPVCoordinates().getVelocity().getX();
                            final double vyn = numState.getPVCoordinates().getVelocity().getY();
                            final double vzn = numState.getPVCoordinates().getVelocity().getZ();
                            final double pxk = kepState.getPVCoordinates().getPosition().getX();
                            final double pyk = kepState.getPVCoordinates().getPosition().getY();
                            final double pzk = kepState.getPVCoordinates().getPosition().getZ();
                            final double vxk = kepState.getPVCoordinates().getVelocity().getX();
                            final double vyk = kepState.getPVCoordinates().getVelocity().getY();
                            final double vzk = kepState.getPVCoordinates().getVelocity().getZ();
                            stream.println( // numState.getDate() +
                                           " " + pxn + " " + pyn + " " + pzn +
                                           " " + vxn + " " + vyn + " " + vzn +
                                           " " + pxk + " " + pyk + " " + pzk +
                                           " " + vxk + " " + vyk + " " + vzk);
                            break;
                            //
                        }
                    }
                }
            }
            System.out.println(resultSavedAs + output1);

        } catch (OrekitException oe) {
            System.err.println(oe.getLocalizedMessage());
            System.exit(1);
        } catch (IOException ioe) {
            System.err.println(ioe.getLocalizedMessage());
            System.exit(1);
        }
    }

    // Step Handler
    private static class StatesHandler implements OrekitFixedStepHandler {

        private final List<SpacecraftState> states;

        /** Simple constructor.
         */
        StatesHandler() {
            // prepare an empty list of SpacecraftState
            states = new ArrayList<SpacecraftState>();
        }

        /** {@inheritDoc} */
        public void handleStep(final SpacecraftState currentState) {

            // add the current state
            states.add(currentState);

        }

        /** Get the list of handled orbital elements.
         * @return orbital elements list
         */
        public List<SpacecraftState> getStates() {
            return states;
        }

    }

}
