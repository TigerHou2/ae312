close all
clear;clc

M1 = 3;
gamma = 1.4;
M4 = 0.5;
theta = deg2rad(10);

p1 = 30;

p0_1 = p1 * (1+(gamma-1)/2*M1^2)^(gamma/(gamma-1));

% ========== Double Compression Corner ==========

% Bend #1

beta1 = TBM_get_beta(M1,theta,gamma);
pres1 = TBM_get_pres(M1,beta1,gamma);

M1n = M1 * sin(beta1);
M2n = sqrt( (M1n^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*M1n^2-1) );
M2 = M2n / sin(beta1-theta);

% Bend #2

beta2 = TBM_get_beta(M2,theta,gamma);
pres2 = TBM_get_pres(M2,beta2,gamma);

M2n = M2 * sin(beta2);
M3n = sqrt( (M2n^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*M2n^2-1) );
M3 = M3n / sin(beta2-theta);

% Normal Shock

M4_comp = sqrt( (M3^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*M3^2-1) );
pres3 = 1 + 2*gamma/(gamma+1)*(M3^2-1);

p4_comp = p1 * pres1*pres2*pres3;

p0_4_comp = p4_comp * (1+(gamma-1)/2*M4_comp^2)^(gamma/(gamma-1));
eff_comp = p0_4_comp / p0_1;

disp(['Double compression corner efficiency:' newline ...
      '    ' num2str(eff_comp)])

% ========== Single Normal Shock ==========
M4_norm = sqrt (( M1^2 + 2/(gamma-1) ) / ( 2*gamma/(gamma-1)*M1^2 - 1 ));
pres4 = 1 + 2*gamma/(gamma+1)*(M1^2-1);

p4_norm = p1 * pres4;

p0_4_norm = p4_norm * (1+(gamma-1)/2*M4_norm^2)^(gamma/(gamma-1));
eff_norm = p0_4_norm / p0_1;

disp(['Normal shock inlet efficiency:' newline ...
      '    ' num2str(eff_norm)])
disp('The double compression corner inlet is more efficient.')