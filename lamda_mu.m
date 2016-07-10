function [DF_DSIG1,MDF_DSIG1,DEP1,DMu_Deffp1,Kp1,TERM1,TERM2,DP,dlamda1,dmu1] = lamda_mu(P1,Q1,S1,DE1,Mu1,STRAN_P_d1,subdEps1)
    
DF_DSIG1 = zeros(6,1);
DP_DSIG = zeros(6,1);
DQ_DSIG = zeros(6,1);
DP = zeros(6,6);
DEP1 = zeros(6,6);
TERM1 = zeros(6,6);
    	
 %    Calculate total stiffness matrix	
	%    Calculate DF_DP
	
	DF_DP = -Mu1; 
	
	%      Calculate DF_DQ
	DF_DQ = 1.0/sqrt(3);
	
	
	%      Calculate DP_DSIG	  
	for K= 1:3
	DP_DSIG(K,1) = 1.0/3.0;
	end
	
	%      Calculate DQ_DSIG	  
	for K= 1:3
	DQ_DSIG(K,1) = (3.0/(2.0*Q1))*S1(K,1);
	end
	for K= 4:6
	DQ_DSIG(K,1) = (3.0/(2.0*Q1))*2.0*S1(K,1);
	end
	
	%    CALCULATE DF_DSIG
	
	for K= 1:6
	DF_DSIG1(K,1) = DF_DP*DP_DSIG(K,1)+DF_DQ*DQ_DSIG(K,1);           
	end
 	
       %    CALCULATE norm of DF_DSIG

        MDF_DSIG1 = sqrt(DF_DSIG1(1,1)^2+DF_DSIG1(2,1)^2+DF_DSIG1(3,1)^2+DF_DSIG1(4,1)^2+DF_DSIG1(5,1)^2+DF_DSIG1(6,1)^2);
        
        %    CALCULATE trace of DF_DSIG
        trDF_DSIG1 = DF_DSIG1(1,1)^2+DF_DSIG1(2,1)^2+DF_DSIG1(3,1)^2+DF_DSIG1(4,1)^2+DF_DSIG1(5,1)^2+DF_DSIG1(6,1)^2;
        
        %Calculate the deviotric part of DF_DSIG
        DF_DSIG_v1 = DF_DSIG1(1,1)+DF_DSIG1(2,1)+DF_DSIG1(3,1);
        DF_DSIG_d1 = sqrt(3/2)*sqrt((DF_DSIG1(1,1)-DF_DSIG_v1/3)^2+(DF_DSIG1(2,1)-DF_DSIG_v1/3)^2+(DF_DSIG1(3,1)-DF_DSIG_v1/3)^2+(DF_DSIG1(4,1))^2+(DF_DSIG1(5,1))^2+(DF_DSIG1(6,1))^2);
	
	%    CALCULATE CONTRIBUTION OF HARDENING
    c1 = 0.7;
    c2 = -10;
    c3 = 0.4;
    DMu_Deffp1 = -c1*c2*c3*(STRAN_P_d1+0.0001)^(c3-1)*exp(c2*(STRAN_P_d1)^c3);
	Kp1 = trDF_DSIG1*DMu_Deffp1*P1;
    
	%   Calculate plastic stiffness matrix
	TERM1 = DE1*(DF_DSIG1)*DF_DSIG1'*DE1;
	TERM2 = DF_DSIG1'*DE1*DF_DSIG1;	
 
	for I=1:6
	for J=1:6
	DP(I,J)=TERM1(I,J)/(Kp1+TERM2);
	end
    end
    
    dlamda1 = DF_DSIG1'*DE1*subdEps1/(DF_DSIG1'*DE1*DF_DSIG1 + Kp1);
 
    % Calculation of dmu
     
    B1 = (1/sqrt(2))*DMu_Deffp1*MDF_DSIG1;
    dmu1 = dlamda1*B1;
	
	
end