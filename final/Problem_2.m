%% Problem 2

close all
clear;clc

disp('Problem 2')
disp(' ')

atm = 101.3e3; % kPa
p0 = 3 * atm;
pb = 1 * atm;
gamma = 1.4;
R = 287;

Aea = 1e-4; % m^2

Atb = 1e-4; % m^2
Aeb = 3e-4; % m^2

gamm = (gamma-1)/2;
gamp = (gamma+1)/2;

Tb = 23 + 273;
T0 = Tb;

% ------------------------------
% Part a

% calculate mdot for nozzle A
% --- first, we find that the flow is choked
ps = p0 * (1+gamm)^(-gamma/2/gamm);
MeA = 1;
% --- 
pe = p0 * (1+(gamma-1)/2*MeA^2)^(-gamma/(gamma-1));
Te = T0 * (1+(gamma-1)/2*MeA^2)^(-1);
ve = MeA * sqrt(gamma*R*Te);
rho_e = pe/R/Te;
mdotA = rho_e * Aea * ve;

% calculate mdot for nozzle B
fun = @(M) 1/M^2 * (1/gamp * (1+gamm * M^2))^(gamp/gamm) - (Aeb/Atb)^2;
Me3 = fzero(fun, [1e-16,1]);
Me6 = fzero(fun, [1,20]);
Pb3 = (1+gamm*Me3^2)^(-gamma/2/gamm) * p0;
Pb6 = (1+gamm*Me6^2)^(-gamma/2/gamm) * p0;
Pb5 = Pb6 * (1+gamma/gamp*(Me6^2-1));
% --- first, we find that the nozzle is over-expanded
% --- from the cals above. Therefore Me = Me6
MeB = Me6;
mdotB = Aeb*p0/sqrt(T0)*sqrt(gamma/R)*MeB/((1+gamm*MeB^2)^(gamp/gamm/2));

if mdotA>mdotB
    disp('a) Nozzle A delivers larger mass flow rate')
elseif mdotB>mdotA
    disp('a) Nozzle B delivers larger mass flow rate')
else
    disp('a) Both nozzles deliver the same mass flow rate')
end
disp(' ')

% ------------------------------
% Part b
if MeB > MeA
     disp('b) TRUE.  ("The Mach number at the exit of Nozzle B is larger than that for Nozzle A")')
else disp('b) FALSE. ("The Mach number at the exit of Nozzle B is larger than that for Nozzle A")')
end
if MeB > 1
     disp('   TRUE.  ("The flow at the exit of Nozzle B is supersonic")')
else disp('   FALSE. ("The flow at the exit of Nozzle B is supersonic")')
end
if Pb6 == 1*atm
     disp('   TRUE.  ("The pressure at the exit of Nozzl B is 1 atm")')
else disp('   FALSE. ("The pressure at the exit of Nozzl B is 1 atm")')
end
     disp('   TRUE.  ("The flow within both nozzles is isentropic")')
     disp('   FALSE. ("The pressure at the inlet of Nozzle A is higher than 3.0 atm")')
disp(' ')