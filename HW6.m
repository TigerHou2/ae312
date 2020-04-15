%% Problem 1

close all
clear;clc

p1 = 100;
gamma = 1.4;
theta1 = deg2rad(-12);
theta2 = deg2rad(8);
mu1 = deg2rad(20);

% ----- Part (a) -----

M1 = 1/sin(mu1);

% ----- Part (b) -----

syms m2 m3 real

eqn = theta1 == PM(m2,gamma) - PM(M1,gamma);
M2 = double(vpasolve(eqn,m2));

eqn = theta2 == PM(m3,gamma) - PM(M2,gamma);
M3 = double(vpasolve(eqn,m3));