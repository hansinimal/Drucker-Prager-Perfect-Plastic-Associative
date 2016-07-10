function [P1,Q1,S1] = PQ(STRESS1)
 %******* Calculate stress invariants p and q
 
    S1 = zeros(4,1);
    
    P1 = (STRESS1(1,1)+STRESS1(2,1)+STRESS1(3,1))/3.0 ;
    
    % Calculate effective stress
    
    S1(1,1) = STRESS1(1,1)-P1;
    S1(2,1) = STRESS1(2,1)-P1;
    S1(3,1) = STRESS1(3,1)-P1;
    S1(4,1) = STRESS1(4,1);
    
    Q1 = sqrt(1/2.0)*sqrt(S1(1,1)^2+S1(2,1)^2+S1(3,1)^2+2*S1(4,1)^2);
    
    
 end
 