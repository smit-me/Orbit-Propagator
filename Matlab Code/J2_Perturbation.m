function [ OUT1 ] = J2_Perturbation( P )
%To find the Gravitational Acceleration on the satellite at Position P. 
%The Earth is assumed to be centred at the Origin.

%OUT1 = Acceleration due to Gravity

%     G = 6.67428 * ( 10 ^ ( -11 ) ) ;
%     Mass_Earth = 5.972 * ( 10 ^ ( 24 ) ) ;
    R_e = 6378136.3 ;
%     J2 = 1.75553e25 ;
%     J2 = -4.84165390000000e-04;
    J2 = 1.082635854e-3;
    GM = 398600.4418e+9;
    
    X = P( 1 ) ;
    Y = P( 2 ) ;
    Z = P( 3 ) ;
    R = norm( P ) ;
    FG = ( GM * R_e * R_e * J2 ) ;
    FGnew = FG * ( 1.5 / ( R ^ 5 ) ) ; 
    F( 1 , 1 ) = ( 5 * X * ( ( Z / R ) ^ 2 ) ) - ( X ) ;
    F( 2 , 1 ) = ( 5 * Y * ( ( Z / R ) ^ 2 ) ) - ( Y ) ;
    F( 3 , 1 ) = ( 5 * Z * ( ( Z / R ) ^ 2 ) ) - ( 3 * Z ) ;
    
%     j2acc = (-1.5) * J2 * (GM / (R^2)) * ((R_e/R)^2) * [ (1 - (5 * ((Z/R)^2))) * (X/R);
%         (1 - (5 * ((Z/R)^2))) * (Y/R);
%         (3 - (5 * ((Z/R)^2))) * (X/R)];
    OUT1 = ( F * FGnew ) ;
%     OUT1 = j2acc;

end

