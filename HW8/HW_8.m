%% Problem 2
clc; clear; close all

alpha = deg2rad(4); 
t_over_c = 0.05;
M_1 = 0;
M_2 = 2;
M_sub = linspace(M_1, 1, 500);
M_sup = linspace(1, M_2, 500);

syms x_c
y_c = 0.385*t_over_c*(1 - 2*x_c)*sqrt(1 - (2*x_c)^2);
phi = diff(y_c,x_c);
cp_upper = 2*(alpha - phi)/sqrt(M_2^2 - 1);
cp_lower = 2*(alpha + phi)/sqrt(M_2^2 - 1);

figure
hold on; grid on; grid minor;
fplot(cp_upper,'r-','LineWidth', 2)
fplot(cp_lower,'k-','LineWidth', 2)
xlim([-0.5 0.5]); ylim([-1 1]);
xlabel('x/c')
ylabel('c$_p$','fontsize', 22)
legend('Upper Surface','Lower Surface')

Cl_0 = 2*pi*(1 + 0.77*t_over_c)*sin(alpha);
Cl_sub = Cl_0./sqrt(1 - M_sub.^2);
Cl_sup = 4*alpha./sqrt(M_sup.^2 - 1);

figure
hold on; grid on; grid minor;
plot(M_sub,Cl_sub,'r-','LineWidth', 2);
plot(M_sup,Cl_sup,'k-','LineWidth', 2);
xlabel('Mach Number')
ylabel('c$_l$','fontsize', 22)
legend('Subsonic','Supersonic')
latexify

%% Functions
function val = nu(M)
gamma = 1.4;
val = sqrt((gamma+1)/(gamma-1))*atan2(sqrt((gamma-1)*(M.^2-1)),sqrt(gamma+1)) - atan2(sqrt(M.^2-1),1);
end