close all
clear;clc

M = 2;
gamma = 1.4;
T = 277.8;
theta = deg2rad(10);

beta = TBM_get_beta(M,theta,gamma);

sin_ang = beta-theta;

hypo = 1/sin(beta);
dist = 1/tan(beta);

d2_d1 = hypo * sin(sin_ang);

disp(['The lower bend is ' num2str(dist) '*D1 to the right of the top bend.'])
disp(['D2 = ' num2str(d2_d1) '*D1'])