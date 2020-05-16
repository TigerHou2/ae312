%% solve for values symbolically

close all
clear;clc

syms M1 T1 p1 gamma positive

T2 = T1 * (1+2*gamma/(gamma+1)*(M1^2-1)) ...
        * ((2+(gamma-1)*M1^2)/((gamma+1)*M1^2));
p2 = p1 * (1+2*gamma/(gamma+1)*(M1^2-1));
M2 = sqrt((M1^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1^2-1));

