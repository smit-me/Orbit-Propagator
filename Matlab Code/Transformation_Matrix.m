%% Input Arguments
    % date = Julian Date of the desired ECI frame
    
%% Output Arguments
    % M = Transfromation matrix that converts J2000 ECI frame to ECI
    %     frame at date time.

function M = Transformation_Matrix(date)

    T = (date - 2451545) / 36525;
    
    %% Tranformation Matrix due to Axial Precession
    
    zeta = 2306.2181 * T + 0.30188 * T^2 + 0.017998 * T^3;
    zeta = deg2rad(zeta / 3600);
    z = 2306.2181 * T + 1.09468 * T^2 + 0.018203 * T^3;
    z = deg2rad(z / 3600);
    theta = 2004.3109 * T - 0.42655 * T^2 - 0.041833 * T^3;
    theta = deg2rad(theta / 3600);
    D = [cos(z) * cos(theta) * cos(zeta) - sin(z) * sin(zeta), -cos(z) * cos(theta) * sin(zeta) - sin(z) * cos(zeta), -cos(z) * sin(theta);
        sin(z) * cos(theta) * cos(zeta) + cos(z) * sin(zeta), -sin(z) * cos(theta) * sin(zeta) + cos(z) * cos(zeta), -sin(z) * sin(theta);
        sin(theta) * cos(zeta), -sin(theta) * sin(zeta), cos(theta)];
    
    %% Transformation Matrix due to Earth's Nutation
    
    [angles, rates] = earthNutation(date); 
    delta_psi = angles(1);
    delta_epsilon = angles(2);
    epsilon_bar = 84381.448 - 46.8150 * T - 0.00059 *T^2 + 0.001813 * T^3;
    epislon_bar = deg2rad(epsilon_bar / 3600);
    epsilon = epsilon_bar + delta_epsilon;
    C = [cos(delta_psi), -sin(delta_psi) * cos(epsilon_bar), -sin(delta_psi) * sin(epsilon_bar);
        cos(epsilon) * sin(delta_psi), cos(epsilon) * cos(delta_psi) * cos(epsilon_bar) + sin(epsilon) * sin(epsilon_bar), cos(epsilon) * cos(delta_psi) * sin(epsilon_bar) - sin(epsilon) * cos(epsilon_bar);
        sin(epsilon) * sin(delta_psi), sin(epsilon) * cos(delta_psi) * cos(epsilon_bar) - cos(epsilon) * sin(epsilon_bar), sin(epsilon) * cos(delta_psi) * sin(epsilon_bar) + cos(epsilon) * cos(epsilon_bar)];
    
    %% Transformation Matrix due to Polar Motion
    
    angularDisp = polarMotion(date);
    xp = angularDisp(1);
    yp = angularDisp(2);
    A = [1, 0, xp;
        0, 1, -yp;
        -xp, yp, 1];
    
    %% Combining both transformation matrices
    
    M = A * C * D;
    
end

