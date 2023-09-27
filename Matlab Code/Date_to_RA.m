function [ OUT ] = Date_to_RA( IN ) 
%To Calculate the Initial Right Ascension of the satellite based on the
%Launch date. The input argument is the number of days past since the last
%Vernal Equinox i.e. 21 March

    IN2 = IN - floor( IN ) ; 
    Angle = ( floor( IN ) * 360 ) / ( 365.25 ) + ( ( IN2 ) * 360 ) + ( 180 ) + 82.5 ;
    Angle = ( Angle ) * ( pi / 180 ) ;
    OUT = mod( Angle , ( 2 * pi ) ) ;
end

