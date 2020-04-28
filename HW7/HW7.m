%% Problem 1
close all
clear;clc

disp('Problem 1')

gamma = 1.4;
R = 287;
A = pi*(0.01/2)^2; % m^2, exit area
Me = 1;
pb = 15e3; % kPa, reservoir pressure
po = pb; % kPa, stagnation pressure
Tb = 20 + 273.15; % K, reservoir temperature
To = Tb; % K, stagnation temperature

pamb = 100; % kPa, ambient pressure

% --- part a) ---
    % flow is choked because 0.5283*pb > pamb
    % flow is sonic at throat
pe = po * (1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1));
Te = To * (1+(gamma-1)/2*Me^2)^(-1);
ve = Me * sqrt(gamma*R*Te);
disp(['a)   Pe = ' num2str(pe) ' kPa' newline ...
      '     Ve = ' num2str(ve) ' m/s'])
  
% --- part b) ---
rho_e = pe/R/Te*1000; % includes conversion from kPa to Pa
mdot = rho_e * A * ve;
disp(['b)   mdot = ' num2str(mdot) ' kg/s'])

% --- part c) ---
F = mdot * ve + (pb-pe) * A * 1000;
disp(['c)   F = ' num2str(F) ' N'])

% --- part d) ---


%% Problem 2
close all
clear;clc

disp('Problem 2')

gamma = 1.4;
R = 287;

% --- part a) ---

% --- part b) ---
M = linspace(0,5,2001);
LHS = sqrt(gamma/R).*M./((1+(gamma-1)/2.*M.^2)).^((gamma+1)/2/(gamma-1));
plot(M,LHS)
grid(gca,'minor')
grid on
xlabel('Mach Number')
ylabel('$\dot{m}\sqrt{T_0}/(p_0A)$')
latexify

% --- part c) ---
LHS_max = max(LHS);
M_max = M;
M_max(LHS~=LHS_max) = [];
disp(['c)   Val = ' num2str(LHS_max) newline ...
      '       M = ' num2str(M_max)])

%% Problem 3
close all
clear;clc

gamma = 1.4;
R = 287;
T0 = 310; % K
p0 = 810; % kPa
pe = 101.3; % kPa
mdot = 1; % kg/s

% --- part a) ---
Mt = 1; % Mach number at throat
A = ( sqrt(gamma/R).*Mt./((1+(gamma-1)/2.*Mt.^2)).^((gamma+1)/2/(gamma-1)) ...
    / mdot / sqrt(T0) * p0 * 1000) ^ (-1);
disp(['a)   A_throat = ' num2str(A) ' m^2'])

% --- part b) ---
syms mm positive
eqn = pe/p0 == (1+(gamma-1)/2*mm^2)^(-gamma/(gamma-1));
Me = double(solve(eqn,mm));
disp(['b)   M_exit   = ' num2str(Me)])

% --- part c) ---
Te = T0 * (1+(gamma-1)/2*Me^2)^(-1);
ve = Me * sqrt(gamma*R*Te);
disp(['c)   V_exit   = ' num2str(ve) ' m/s'])

%% Problem 4
close all
clear;clc

gamma = 1.3;
cp = 1.2;

% stand-in variables, won't affect results
R  = 276.923;
At = 1;
Ae = 1.5 * At;

T0 = 1670; % K
p0 = 3e3; % kPa
pb = 101.3;

% --- converging nozzle case ---
    % flow is choked because exit pressure lower than 0.5283 * p0
    % flow is Mach 1 at exit (throat)
Mt = 1;
pt = p0 * (1+(gamma-1)/2*Mt^2)^(-gamma/(gamma-1));
Tt = T0 * (1+(gamma-1)/2*Mt^2)^(-1);
vt = sqrt(gamma*R*Tt);
rho_t = pt/R/Tt * 1000; % includes conversion from kPa to Pa
mdot_t = rho_t * At * vt;
F_orig = mdot_t * vt + (pt-pb) * At * 1000;

% --- cd nozzle case ---
options = optimoptions('fmincon','Display','off');
Me = fmincon(@(mm) abs(mdot_t/Ae/p0/1000*sqrt(T0) - sqrt(gamma/R)*mm...
         /(1+(gamma-1)/2*mm^2)^((gamma+1)/2/(gamma-1))), 2, -1, -1.5, ...
         [], [], [], [], [], options);
% syms mm positive
% eqn = pe/p0 == (1+(gamma-1)/2*mm^2)^(-gamma/(gamma-1));
% Me = double(solve(eqn,mm));
Te = T0 * (1+(gamma-1)/2*Me^2)^(-1);
ve = Me * sqrt(gamma*R*Te);
pe = p0 * (1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1));
rho_e = pe/R/Te * 1000; % includes conversion from kPa to Pa
mdot_e = rho_e * Ae * ve;
F_cd = mdot_e * ve + (pe-pb) * Ae * 1000;

disp(['Thrust increase = ' num2str((F_cd-F_orig)/F_orig*100) '%'])
