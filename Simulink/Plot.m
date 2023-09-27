load_data = load('elts_pv.txt');
Orekit_data = load_data';
pos_orekit = Orekit_data(1:3,:);
vel_orekit = Orekit_data(4:6,:);
poskep_orekit = Orekit_data(7:9,:);
velkep_orekit = Orekit_data(10:12,:);
ttpos= timeseries2timetable(Sim_pos_exp);
ttvel= timeseries2timetable(Sim_vel_exp);
pos_sim = ttpos.Data';
vel_sim = ttvel.Data';


%calculate relative position 
OrekitRelativePosition = vecnorm(pos_mat - pos_orekit,2,1); % matlab wrt orekit 
SimulinkRelativePosition = vecnorm(pos_sim - pos_orekit,2,1); % sim wrt orekit 

plottime = 0:60:(Final_Time+60);
%plot relative position wrt Keplerian 
f1 = figure ;
OrekitRelativePositionKm = OrekitRelativePosition/1000;
SimulinkRelativePositionKm = SimulinkRelativePosition/1000;
plot(plottime,OrekitRelativePositionKm,plottime,SimulinkRelativePositionKm)
title('Positional Error Comparison')
xlabel("Time")
ylabel("Relative position (km)")
legend("SGP4","Spherical Harmonics")
xlim([0 172800])
xticks(0:43200:172800)
xticklabels({'22/8 0:00','22/8 12:00','23/8 0:00', '23/8 12:00', '24/8 0:00'})


%calculate relative velocity 
OrekitRelativeVelocity = vecnorm(vel_mat - vel_orekit,2,1);
SimulinkRelativeVelocity = vecnorm(vel_sim - vel_orekit,2,1);

%plot relative velocity wrt Keplerian
f2 = figure ;
plot(plottime,OrekitRelativeVelocity,plottime,SimulinkRelativeVelocity)
title('Velocity Error Comparison')
xlabel("Time")
ylabel("Velocity deviation (m/s)")
legend("SGP4","Spherical Harmonics")
xlim([0 172800])
xticks(0:43200:172800)
xticklabels({'22/8 0:00','22/8 12:00','23/8 0:00', '23/8 12:00', '24/8 0:00'})


f3 = figure ;
plot(plottime,(pos_mat - pos_orekit)/1000)
title('SGP4 vs. Two-body Keplerian')
xlabel("Time")
ylabel("Position deviation (km)")
legend("\DeltaX_{ICRF}","\DeltaY_{ICRF}","\DeltaZ_{ICRF}" )
xlim([0 172800])
xticks(0:43200:172800)
xticklabels({'22/8 0:00','22/8 12:00','23/8 0:00', '23/8 12:00', '24/8 0:00'})



f4 = figure ;
plot(plottime,(pos_sim - pos_orekit)/1000)
title('Spherical Harmonics vs.Two Body Keplerian')
xlabel("Time")
ylabel("Position deviation (km)")
legend("\DeltaX_{ICRF}","\DeltaY_{ICRF}","\DeltaZ_{ICRF}" )
xlim([0 172800])
xticks(0:43200:172800)
xticklabels({'22/8 0:00','22/8 12:00','23/8 0:00', '23/8 12:00', '24/8 0:00'})
