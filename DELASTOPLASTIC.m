function [DF_DSIG1,M_EF_DF_DSIG1,DEP1,TERM1,TERM2,DP] = DELASTOPLASTIC(P1,Q1,S1,DE1,Mu1)
    
DF_DSIG1 = zeros(4,1);
DP_DSIG = zeros(4,1);
DQ_DSIG = zeros(4,1);
DP = zeros(4,4);
DEP1 = zeros(4,4);
TERM1 = zeros(4,4);
EF_DF_DSIG1= zeros(4,1);
    	
    %    Calculate total stiffness matrix	
	%    Calculate DF_DP
	
	DF_DP = -Mu1; 
	
	%      Calculate DF_DQ
	DF_DQ = 1;
		
	%      Calculate DP_DSIG	  
	for K= 1:3
	DP_DSIG(K,1) = 1.0/3.0;
	end
	
	%      Calculate DQ_DSIG	  
	for K= 1:3
	DQ_DSIG(K,1) = (1/(2.0*Q1))*S1(K,1);
    end
    
    DQ_DSIG(4,1) = (1/(2.0*Q1))*2*S1(4,1);
	
	%    CALCULATE DF_DSIG
	
	for K = 1:4
	DF_DSIG1(K,1) = DF_DP*DP_DSIG(K,1)+DF_DQ*DQ_DSIG(K,1);           
    end
    
    %    CALCULATE effective DF_DSIG
	
	for K = 1:3
	EF_DF_DSIG1(K,1) = DF_DSIG1(K,1)-(DF_DSIG1(1,1)+DF_DSIG1(2,1)+DF_DSIG1(3,1))/3.0;           
    end
    
	EF_DF_DSIG1(4,1) = DF_DSIG1(4,1);           
		
    %    CALCULATE norm of effective DF_DSIG

    M_EF_DF_DSIG1 = sqrt(EF_DF_DSIG1(1,1)^2+EF_DF_DSIG1(2,1)^2+EF_DF_DSIG1(3,1)^2+EF_DF_DSIG1(4,1)^2);
  

    
	%   Calculate plastic stiffness matrix
	TERM1 = DE1*(DF_DSIG1)*DF_DSIG1'*DE1;
	TERM2 = DF_DSIG1'*DE1*DF_DSIG1;

	for I=1:4
	for J=1:4
	DP(I,J)=TERM1(I,J)/(TERM2);
	end
	end
		
	%   Calculate elasto plastic stiffness matrix	
	
	for I=1:4
	for J=1:4
	DEP1(I,J)=DE1(I,J)-DP(I,J);
	end
	end
	
end