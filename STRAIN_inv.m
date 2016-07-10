function [STRAN_v,STRAN_d] = STRAIN_inv(STRAN)
    % Calculate volumetric strain increment
        STRAN_v = (STRAN(1,1)+ STRAN(2,1) + STRAN(3,1));
     
    % Calculate deviatoric strain increment
	STRAN_d = sqrt(1/2)*sqrt((STRAN(1,1)-STRAN_v/3)^2+(STRAN(2,1)-STRAN_v/3)^2+(STRAN(3,1)-STRAN_v/3)^2+2*(STRAN(4,1))^2 );
end