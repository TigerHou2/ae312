%% solve for values symbolically

close all
clear;clc

syms Us a1 Ms p2 T2 Ur gamma positive

% find M1 in terms of M2 in NS equation
% syms M1 M2 real
% NS = M2^2 == (M1^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*M1^2-1);
% isolate(NS,M1)

% find upstream conditions in incident shock
M2 = Us/a1;
M1 = sqrt((M2^2*(gamma-1)+2)/(2*gamma*M2^2-gamma+1));
T1 = T2 / (1+(gamma-1)/2*M1^2) * (1+(gamma-1)/2*M2^2);
Ug = Us * (1-M2/M1*sqrt(T2/T1));

% find upstream conditions in reflected shock
Mr = Ur/a1;
Mo = sqrt((Mr^2*(gamma-1)+2)/(2*gamma*Mr^2-gamma+1));
Tr = T1;
To = Tr / (1+(gamma-1)/2*Mo^2) * (1+(gamma-1)/2*Mr^2);
eqn = 1 + Ug/Ur == Mo/Mr*sqrt(To/Tr);
isolate(eqn,Ur)