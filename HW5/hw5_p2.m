close all
clear;clc

theta = deg2rad(10);
M = 3;
gamma = 1.4;

beta = TBM_get_beta(M,theta,gamma);
p_ratio_1 = TBM_get_pres(M,beta,gamma);

M1n = M * sin(beta);
M2n = sqrt( (M1n^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*M1n^2-1) );

M2 = M2n / sin(beta-theta);

beta2 = TBM_get_beta(M2,theta,gamma);
p_ratio_2 = TBM_get_pres(M,beta2,gamma);

p_ratio_tot = p_ratio_1 * p_ratio_2;

figure(2)
hold on

x_axis = linspace(0,1,1000);

y_wall = ones(size(x_axis));
y_wall(x_axis> 0.7) = p_ratio_tot;

y_inch = ones(size(x_axis));
y_inch(x_axis>0.5) = p_ratio_1;
y_inch(x_axis>0.8) = p_ratio_tot;

plot(x_axis,y_wall,'LineWidth',2)
plot(x_axis,y_inch,'--','LineWidth',2)

text(  0.05, 1.2, '$p/p_{\infty} = 1$', 'FontSize',16)
text(0.5, p_ratio_1+0.2, ...
       strcat('$p/p_{\infty} = ',num2str(p_ratio_1),'$'), 'FontSize',16)
text(0.7, p_ratio_tot+0.2, ...
       strcat('$p/p_{\infty} = ',num2str(p_ratio_tot),'$'), 'FontSize',16)

hold off

xlabel('Normalized Distance from Wedge Tip','FontSize',16)
ylabel('$p/p_{\infty}$','FontSize',16)

set(gca,'FontSize',16)
grid(gca,'minor')
grid on

legend('Along the Wall','1-in. Away From the Wall',...
       'FontSize',16,'Location','northwest')

latexify