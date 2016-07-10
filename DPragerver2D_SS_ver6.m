% ****************** Drucker Prager 2D ******************%
% ****************** With hardening *********************%

%  Single element model
%  Drained Simple Shear

% Modified Euler method; without substepping 
% Jaccorbian calculated based on the trial stress 
% yield check on trial stress

%%%%%%%%%%%%%%%%Modifications%%%%%%%%%%%%%%%%%%%%%%%%
% stress and strain inavariants
% update of cumulative deviatoric plastic strain
% Kp and B based on modulus of effective df_dsig
% Mu can be calculated based on the equation(c2 = -200) or kept constant

% Associative and perfect plasticity

iteration = 1000;
P_initial = 100;
dEps12 = 0.0001;

K0 = 0.4;

% Storing output variables
P = zeros(iteration,1);
Q = zeros(iteration,1);


P_trial = zeros(iteration,1);
Q_trial = zeros(iteration,1);
YIELD_trial = zeros(iteration,1);
YIELD_initial = zeros(iteration,1);
dlamda_t = zeros(iteration,1);
Kp_t = zeros(iteration,1);

Mu = zeros(iteration,1);
Epsv = zeros(iteration,1);
Epsd = zeros(iteration,1);
STRAN_P_vol = zeros(iteration,1);
STRAN_P_dev = zeros(iteration,1);
deriv_f_Mu = zeros(iteration,1);

S22 = zeros(iteration,1);
S11 = zeros(iteration,1);
S12 = zeros(iteration,1);
S33 = zeros(iteration,1);
Eps11 = zeros(iteration,1);
Eps12 = zeros(iteration,1);

Kp = zeros(iteration,1);
dlamda = zeros(iteration,1);
dmu = zeros(iteration,1);
shear_ratio = zeros(iteration,1);
YIELD = zeros(iteration,1);
stress_rotation = zeros(iteration,1);
strain_rotation = zeros(iteration,1);


%  Defining zero matrices

S = zeros(4,1);
S_trial = zeros(4,1);
STRESS = zeros(4,1);
STRESS_trial = zeros(4,1);
STRAN = zeros(4,1);
STRAN_P = zeros(4,1);
DSTRAN_P = zeros(4,1);

dEps = zeros(4,1);
subdEps = zeros(4,1);

DF_DSIG = zeros(4,1);
DP_DSIG = zeros(4,1);
DQ_DSIG = zeros(4,1);
DE = zeros(4,4);
DP = zeros(4,4);
DEP = zeros(4,4);
TERM1 = zeros(4,4);

DSTRESS = zeros(4,1);
DSTRESS_P = zeros(4,1);
DSTRESS_E = zeros(4,1);


%*******Initialization**************%

    % Initial index
    a = 1;
    
    % Initial stress

    STRESS(1,1) = P_initial + 0.0001; % avoid dividing by zero
    STRESS(2,1) = K0*P_initial;
    STRESS(3,1) = K0*P_initial;

  % Calculate initial parameters
      
	[P(1,1),Q(1,1),S] = PQ(STRESS);
    [G,BULK,DE] = DELASTIC(P(1,1));
    
    STRAN_P = 0;
    STRAN_P_d = 0;
    Mu(1,1) = HARD(STRAN_P_d);
      
    YIELD(1,1) = Q(1,1)-Mu(1,1)*P(1,1);
   
    % Calculation of elastic predictor
       
        dEps =  [  -DE(1,4)*dEps12/DE(1,1);
	               0;
	               0;
	               dEps12];
    
    while a < iteration      
       
      %  Calculate the yield function based on the initial stress
       
	   YIELD_initial(a,1) = Q(a,1)-Mu(a,1)*P(a,1);   
           
       % Calculate the trial stress
       STRESS_trial = STRESS + DE*dEps;
       [P_trial(a,1),Q_trial(a,1),S_trial] = PQ(STRESS_trial);
       %  Calculate the yield function based on the trial stress
 	   YIELD_trial(a,1) = Q_trial(a,1)-Mu(a,1)*P_trial(a,1);   
	
	  if YIELD_trial(a,1) <= -0.00001
        
       % Calculation of Jacoorbian
	   DEP = DE ;
       
      %  Update stress and strain  
	   STRESS = STRESS_trial;
	   STRAN = STRAN + dEps;
       STRAN_P_d = 0;
       STRAN_P_v = 0;
       STRAN_P_dev(a,1) = STRAN_P_d;
    % During elastic range DE(1,4)=0,dEps1=0,dsig1=0 ---only shear stress increase during elastic, Mu and P constant, only Q increase     
	  end 
	
	if YIELD_trial(a,1) > -0.00001
    
    % Calculation of Jacoorbian based on trial stress
    [DF_DSIG_t,M_EF_DF_DSIG_t,DEP,TERM1_t,TERM2_t,DP_t] = DELASTOPLASTIC(P_trial(a,1),Q_trial(a,1),S_trial,DE,Mu(a,1));
   
     % Calculation of dlamda based on trial stress
    dlamda_t(a,1) = YIELD_trial(a,1)/(DF_DSIG_t'*DE*DF_DSIG_t );
    check1(a,1) = DF_DSIG_t'*DE*dEps;
    check2(a,1) = DF_DSIG_t'*DE*DF_DSIG_t;
    
     % Calculation of stress increment
     
    DSTRESS_P_t = dlamda_t(a,1)*DE*DF_DSIG_t;
        
     %*********** Update stress, strain and other variables *******%
        
	%  Update stress and strain
    
	STRESS = STRESS_trial - DSTRESS_P_t;
	STRAN = STRAN + dEps;
    
     % Calculate incremental plastic strain
    DSTRAN_P =  dlamda_t(a,1)*DF_DSIG_t;
    
    % Calculate deviatoric incremental plastic strain 
 
    [DSTRAN_P_v,DSTRAN_P_d] = STRAIN_inv(DSTRAN_P);
    STRAN_P_d = STRAN_P_d + DSTRAN_P_d;
    STRAN_P_v = STRAN_P_v + DSTRAN_P_v;
    STRAN_P_dev(a,1) = STRAN_P_d;
    STRAN_P_vol(a,1) = STRAN_P_v;
      end
       
     a = a +1;
   
    Eps12(a,1) = Eps12(a-1,1) + dEps12;
	   
	[P(a,1),Q(a,1),S] = PQ(STRESS);
   % Calculate elastic stiffness matrix
	[G,BULK,DE] = DELASTIC(P(a,1));
    Mu(a,1) = HARD(STRAN_P_d);
    
    dEps =  [  -DEP(1,4)*dEps12/DEP(1,1);
	0;
	0;
	dEps12];
      
       % storing stress result
  
          S11(a,1) = STRESS(1,1);
          S22(a,1) = STRESS(2,1);
          S33(a,1) = STRESS(3,1);
          S12(a,1) = STRESS(4,1);
          
          Eps11(a,1) = STRAN(1,1);
          Eps22(a,1) = STRAN(2,1);
          Eps33(a,1) = STRAN(3,1);
          Eps12(a,1) = STRAN(4,1);
   
          
          shear_ratio(a,1) = S12(a,1)/S11(a,1);
          stress_rotation(a,1) = 0.5*atan(2*S12(a,1)/(S22(a,1)-S11(a,1)))*180/3.14;
          strain_rotation(a,1) = 0.5*atan(2*Eps12(a,1)/(Eps22(a,1)-Eps11(a,1)))*180/3.14; % for coaxial stress and strain rotation should coincide
	end
  
 
 % plot results
 
 figure1 = figure('Name', '2D Simple Shear- stress and strains')

 subplot(2,3,1,'parent',figure1)
 plot(Eps12,S11)
 title('Eps12 vs S11')
 xlabel('Eps12') % x-axis label
 ylabel('S11') % y-axis label
 
 subplot(2,3,2,'parent',figure1)
 plot(Eps12,S33)
 title('Eps12 vs S33')
 xlabel('Eps12') % x-axis label
 ylabel('S33') % y-axis label
 
 subplot(2,3,3,'parent',figure1)
 plot(Eps12,S22)
 title('Eps12 vs S22')
 xlabel('Eps12') % x-axis label
 ylabel('S22') % y-axis label
 
 subplot(2,3,4,'parent',figure1)
 plot(Eps12,S12)
 title('Eps12 vs S12')
 xlabel('Eps12') % x-axis label
 ylabel('S12') % y-axis label
 
 subplot(2,3,5,'parent',figure1)
 plot(Eps12,Eps11)
 title('Eps12 vs Eff11')
 xlabel('Eps12') % x-axis label
 ylabel('Eff11') % y-axis label
 
 subplot(2,3,6,'parent',figure1)
 plot(Eps12,Eps22)
 title('Eps12 vs Eff22')
 xlabel('Eps22') % x-axis label
 ylabel('Eff22') % y-axis label
 
 figure2 = figure('Name', '2D Simple Shear and Hardening Parameters')
 
 subplot(2,3,1,'parent',figure2)
 plot(Eps12,Mu)
 title('Eps12 vs Mu')
 xlabel('Eps12') % x-axis label
 ylabel('Mu') % y-axis label
 

 subplot(2,3,2,'parent',figure2)
 plot(Eps12,dlamda_t)
 title('Eps12 vs dlamda')
 xlabel('Eps12') % x-axis label
 ylabel('dlamda') % y-axis label

 
 
 figure3 = figure('Name', '2D Simple Shear')
 
 subplot(2,3,1,'parent',figure3)
 plot(Eps12,shear_ratio)
 title('Eps12 vs shear ratio')
 xlabel('Eps12') % x-axis label
 ylabel('shear ratio') % y-axis label
 
 subplot(2,3,2,'parent',figure3)
 plot(Eps12,YIELD_initial)
 title('Eps12 vs YIELD initial')
 xlabel('Eps12') % x-axis label
 ylabel('YIELD initial') % y-axis label
 
	
 subplot(2,3,5,'parent',figure3)
 plot(Eps12,stress_rotation)
 title('Eps12 vs strain rotation')
 xlabel('Eps12') % x-axis label
 ylabel('stress rotation') % y-axis label
 
 subplot(2,3,6,'parent',figure3)
 plot(Eps12,strain_rotation,'color','r')
 hold on;
 plot(Eps12,stress_rotation,'color','b')
 
 title('Eps12 vs rotation')
 legend('strain rotation','stress rotation')
 xlabel('Eps12') % x-axis label
 ylabel('rotation') % y-axis label
 
 figure4 = figure('Name', '2D Simple Shear- stress and strain invariants')
 
 subplot(2,3,1,'parent',figure4)
 plot(P,Q)
 title('P vs Q')
 xlabel('P') % x-axis label
 ylabel('Q') % y-axis label
 
 subplot(2,3,2,'parent',figure4)
 plot(Eps12,Q)
 title('Eps12 vs Q')
 xlabel('Eps12') % x-axis label
 ylabel('Q') % y-axis label
 
 subplot(2,3,3,'parent',figure4)
 plot(Eps12,P)
 title('Eps12 vs P')
 xlabel('Eps12') % x-axis label
 ylabel('P') % y-axis label
 
 subplot(2,3,4,'parent',figure4)
 plot(Eps12,STRAN_P_dev)
 title('Eps12 vs STRAN_P_d')
 xlabel('Eps12') % x-axis label
 ylabel('STRAN_P_d') % y-axis label
 
 subplot(2,3,5,'parent',figure4)
 plot(Eps12,STRAN_P_vol)
 title('Eps12 vs STRAN_P_v')
 xlabel('Eps12') % x-axis label
 ylabel('STRAN_P_v') % y-axis label
	