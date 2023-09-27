function OUT = SGP_acc(IN1,UTC)

% IN1 - Position Vector in ECI
% OUT - Acceleration Vector in ECI 

P = IN1 ;  
Acc_gravity = Gravity(P);
Acc_Oblateness = J2_Perturbation( P ) ;
Acc = Acc_gravity + Acc_Oblateness ;
  
OUT = Acc ;
end