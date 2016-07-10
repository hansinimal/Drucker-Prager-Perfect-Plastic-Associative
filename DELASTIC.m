function [G1,BULK1,DE1] = DELASTIC(P1)
    
        DE1 = zeros(4,4);
    
        G1 = 1.25*10^4;
	BULK1 = 2.5*10^4;
	
	for K1 = 1: 3
	for K2 = 1: 3
	DE1(K2,K1) = BULK1-2.0*G1/3.0;
	end 
	DE1(K1,K1) = BULK1+4.0*G1/3.0;
	end 
	
	DE1(4,4) = G1 ;
	
end