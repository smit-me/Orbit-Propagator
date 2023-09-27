%Initialize simulink model 
mission.StartDate = datetime(2019, 8, 22, 0, 0, 0);
mission.Duration  = hours(48);

%Specify Keplerian orbital elements for the satellite
mission.Satellite.SemiMajorAxis  = 7100000.00; % meters
mission.Satellite.Eccentricity   = 0.05;
mission.Satellite.Inclination    = 30;    % deg
mission.Satellite.RAAN           = 50;    % deg
mission.Satellite.ArgOfPeriapsis = 0;    % deg
mission.Satellite.TrueAnomaly    = 0;   % deg

%create a new satelliteScenario object 
mission.mdl = "OrbitPropagatorBlock";
open_system(mission.mdl);

%Define the path to the Orbit Propagator block in the model.
mission.Satellite.blk = mission.mdl + "/Orbit Propagator";

%Set satellite initial conditions
set_param(mission.Satellite.blk, ...
    "startDate",      num2str(juliandate(mission.StartDate)), ...
    "stateFormatNum", "Orbital elements", ...
    "orbitType",      "Keplerian", ...
    "semiMajorAxis",  "mission.Satellite.SemiMajorAxis", ...
    "eccentricity",   "mission.Satellite.Eccentricity", ...
    "inclination",    "mission.Satellite.Inclination", ...
    "raan",           "mission.Satellite.RAAN", ...
    "argPeriapsis",   "mission.Satellite.ArgOfPeriapsis", ...
    "trueAnomaly",    "mission.Satellite.TrueAnomaly");

%Set the position and velocity output frame
set_param(mission.Satellite.blk, ...
    "centralBody",  "Earth", ...
    "outportFrame", "ICRF");

% Configure the propagator
set_param(mission.Satellite.blk, ...
    "propagator",   "Numerical (high precision)", ...
    "gravityModel", "Spherical Harmonics", ...
    "earthSH",      "EGM2008", ... % Earth spherical harmonic potential model
    "shDegree",     "120", ... % Spherical harmonic model degree and order
    "useEOPs",      "on", ... % Use EOP's in ECI to ECEF transformations
    "eopFile",      "aeroiersdata.mat"); % EOP data file

% Apply model-level solver setting
set_param(mission.mdl, ...
    "SolverType", "Fixed-step", ...
    "SolverName", "FixedStepAuto", ...
    "RelTol",     "1e-6", ...
    "AbsTol",     "1e-7", ...
    "StopTime",   string(seconds(mission.Duration)));

%Save model output port data as a dataset of time series objects.
set_param(mission.mdl, ...
    "SaveOutput", "on", ...
    "OutputSaveName", "yout", ...
    "SaveFormat", "Dataset");

%Simulate and export data for Plot
mission.SimOutput = sim(mission.mdl);
mission.Satellite.TimeseriesPosECEF = mission.SimOutput.yout{1}.Values; %Extract position data
mission.Satellite.TimeseriesVelECEF = mission.SimOutput.yout{2}.Values; %Extract velocity data
mission.Satellite.TimeseriesPosECEF.TimeInfo.StartDate = mission.StartDate; %Initialize start data
mission.Satellite.TimeseriesVelECEF.TimeInfo.StartDate = mission.StartDate;
Sim_pos_exp = mission.Satellite.TimeseriesPosECEF; %Save position data for plot 
Sim_vel_exp = mission.Satellite.TimeseriesVelECEF; %Save veloctiy data for plot 


