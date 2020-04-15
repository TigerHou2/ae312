close all
clear;clc

gamma = 1.4;

figure(1)
hold on

theta = linspace(0,60,120001);
M = [1, 1.2, 1.4, 1.6, 1.8, 2, 2.4, 2.8, 3.2, 4, 5];

bw = zeros(size(theta));
bs = zeros(size(theta));

for i = 1:length(M)
    [bw,bs] = TBM_get_beta(M(i),deg2rad(theta),gamma);
    theta_temp = theta;
    theta_temp(imag(bw)~=0) = [];
    bw(imag(bw)~=0) = [];
    bs(imag(bs)~=0) = [];
    bw = rad2deg(bw);
    bs = rad2deg(bs);
    plot(theta_temp,bw,'Color','k')
    plot(theta_temp,bs,'Color','k')
    if length(bw) >= 2
        text(theta_temp(2),bw(2),strcat('M=',num2str(M(i))),...
             'FontSize',10)
    end
end

hold off

ylim([0,90])
set(gca,'FontSize',16)
grid(gca,'minor')
grid on

xlabel('$\theta$ (${}^o$)','FontSize',16)
ylabel('$\beta$ (${}^o$)' ,'FontSize',16)

latexify