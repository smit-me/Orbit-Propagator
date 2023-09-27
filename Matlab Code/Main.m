%SGP4 - Orbit Propagator
format long ;
disp( datetime( 'now' ) ) ;
tic ;

%General Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Mass_Satellite = 1 ;                                                        %Satellite Mass in Kg
Radius_Earth = 6378000 ;   %6378kms                                         %Earth's radius in metres
Radius_Satellite = 722000 ; %722kms  total= 7100kms                         %Satellite altitude in metres
R_Initial = Radius_Earth + Radius_Satellite ;                               %Distance from Earth's Center to Satellite COM in metres
Ang_Vel_Earth = ( 360 ) / ( 86164 ) ;                                       %Angular Velocity of Earth in Degrees/sec.
Launch_Date = decyear( 2019 , 8 , 22 , 0 , 0 , 0 ) ;                        %All times in IST (Indian Standard Time)
Last_Vernal_Equinox = decyear( 2019 , 3 , 21 , 17 , 30 , 0 ) ;              %Assuming Exact Vernal Equinox at 12 noon at Greenwich Meridian
Launch_Difference = ( Launch_Date - Last_Vernal_Equinox ) * ( 365 ) ;       %Number of days between Launch date and last Vernal Equinox i.e. 21 March

% Ground Station Parameters.
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Latitude_GS = 28.361627 ;
% Longitude_GS = 75.591374 ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial Orbital Elements can come here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RA_Satellite = 50 ; %95.2562   Date_to_RA( Launch_Difference )            % Right Ascension of Satellite in Radians; Depends on the Date.
% Perigee_Argument = ( 0 ) ; % 93.4872                                      % Argument of Periapsis in Radians. Constant for a particular Launch Site.
% Eccentricity = 0.05 ;                                                     % Eccentricity
% Mean_Anomaly = 0 ;                                                        % Mean Anomaly in Radians.
% Semi_Major_Axis = 7100 ;                                                  % Semi Major Axis in Metres.
% Inclination = ( 30 * ( pi / 180 ) ) ;                                     % Inclination in Radians.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Initial_P = [4335.602427335707; 5166.969768837508; 0] * 10^3;               % Initial Position in J2000 frame
Initial_V = [-5.22585263; 4.385011015; 3.93860623] * 10^3;                  % Initial Velocity in J2000 frame
M = Transformation_Matrix(juliandate(datetime('now')));                     % M is the transformation Matrix which converts J2000 ECI frame to given time ECI frame
Positions( : , 1 ) = M * Initial_P ;                                        % Initial Position in current ECI frame
Velocities( : , 1 ) = M * Initial_V ;                                       % Initial Velocity in current ECI frame


%Propagator Conditions and Limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sampt = 60 ;                                                                 %Sampling Time in Seconds
Time = 0 ;
Counter = 0 ;
Final_Time = 172740 ;                                                        %Total Time in Orbit in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = zeros( 1 , ( Final_Time / Sampt ) + 1 ) ;

for Time = 0 : Sampt : Final_Time

    Counter = Counter + 1 ;
    
%To find Position and Velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [ Positions( : , Counter + 1 ) , Velocities( : , Counter + 1 ) ] = SGP_Model_RK4(Positions(:, Counter), Velocities(:, Counter), Sampt);
    
end

pos_mat = inv(M) * Positions;     % Converting Positions to J2000 ECI frame
vel_mat = inv(M) * Velocities;   % Converting Velocities to J2000 ECI frame

disp('End of Calculation');
toc ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
