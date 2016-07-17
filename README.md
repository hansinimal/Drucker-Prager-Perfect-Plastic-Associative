"# DruckerP-model" 
% ****************** Drucker Prager 2D ******************%
% ****************** With hardening *********************%

%  Single element model
%  Drained Simple Shear

% Modified Euler method; without substepping 
% Jaccorbian calculated based on the trial stress 
% yield check on trial stress

%%%%%%%%%%%%%%%% Modifications%%%%%%%%%%%%%%%%%%%%%%%%
% stress and strain inavariants
% update of cumulative deviatoric plastic strain
% Kp and B based on modulus of effective df_dsig
% Mu can be calculated based on the equation(c2 = -200) or kept constant

% Associative and perfect plasticity
