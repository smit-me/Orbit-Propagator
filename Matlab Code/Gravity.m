function [ OUT1 ] = Gravity( P )
%To find the Gravitational Acceleration on the satellite at Position P. 

%OUT1 = Acceleration due to Gravity

%     G = 6.67428 * ( 10 ^ ( -11 ) ) ;
%     Mass_Earth = 5.972 * ( 10 ^ ( 24 ) ) ;
    GM = 398600.4418e+9;
    P_Norm = norm( P ) ;
    FG1 = ( GM ) / ( P_Norm ^ 3 ) ;
    OUT1 = ( -1 ) * ( P * FG1 ) ;

end

