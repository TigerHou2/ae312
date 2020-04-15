close all
clear;clc

M1 = 2.4;
M2 = 1.6;
gamma = 1.4;

theta = deg2rad(15);

alpha = deg2rad(8.3315);

theta1 = theta - alpha;
theta2 = alpha;

beta1 = TBM_get_beta(M1,theta1,gamma);
beta2 = TBM_get_beta(M2,theta2,gamma);

p3 = TBM_get_pres(M1,beta1,gamma);
p4 = TBM_get_pres(M2,beta2,gamma);

disp(['alpha = ' num2str(rad2deg(alpha)) ' degrees (clockwise)'])