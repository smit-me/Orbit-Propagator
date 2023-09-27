function [ OUT1 , OUT2 , OUT3] = SGP_Model_RK4( P , V , dt )

% OUT1 - New position vector
% OUT2 - New velocity vector 

% Assumptions 
% 1 - Disturbance forces due air drag and SRP remain constant during integration
% 2 - Position of the Sun remains constant for the time step
          
    Pnew = P ;
    Vnew = V ;
        

       
    kv1 = SGP_acc(Pnew) ;  
    kp1 = Vnew ; 

    kv2 = SGP_acc(Pnew + kp1*(dt/2)) ;  
    kp2 = Vnew + kv1*(dt/2) ;

    kv3 = SGP_acc(Pnew + kp2*(dt/2)) ;  
    kp3 = Vnew + kv2*(dt/2) ; 

    kv4 = SGP_acc(Pnew + kp3*dt) ;  
    kp4 = Vnew + kv3*dt ; 

    Vnew = Vnew + (1/6)*(dt)*( kv1 + 2*kv2 + 2*kv3 + kv4 ) ; 
    Pnew = Pnew + (1/6)*(dt)*( kp1 + 2*kp2 + 2*kp3 + kp4 ) ; 
        

    
    OUT1 = Pnew ;
    OUT2 = Vnew ;
    OUT3 = kv1;
end


