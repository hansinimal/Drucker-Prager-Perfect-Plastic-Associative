"# DruckerP-model" 
% ****************** Drucker Prager 2D ******************%
% ****************** Without hardening *********************%
% ***********Associative and perfect plasticity ***********%

%  Boundary conditions********

%  Single element model
%  Drained Simple Shear

% Numerical method ***************

% Single step backward Euler method; without substepping 
% Jaccorbian calculated based on the trial stress 
% yield check on trial stress
% DF_DSIG, DEP, dlamda based on trial stress
% dlamda = yield trial/()trial

%%%%%%%%%%%%%%%% Modifications%%%%%%%%%%%%%%%%%%%%%%%%
% stress and strain inavariants-p and q
% update of cumulative deviatoric plastic strain
% Kp and B based on modulus of effective df_dsig
% Mu can be calculated based on the equation(c2 = -200) or kept constant

