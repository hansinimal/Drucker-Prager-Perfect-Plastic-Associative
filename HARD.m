function [Mu1] = HARD(STRAN_P_d1)
    a = 0.7;
    b = -200;
    c = 0.4;
  if STRAN_P_d1 <= 0.0001
   STRAN_P_d1 = 0.0001;
  end
  
    Mu1 = a*(1-exp(b*(STRAN_P_d1)^c));
   % Mu1 = 0.7;
  
end