function [cpu,cpl] = cp_Joukowski(x,alpha,tc,M)
%CP_JOUKOWSKI Summary of this function goes here
%   Detailed explanation goes here

c = max(x)-min(x);
t = c * tc;
theta = c .* ( (77.*t.*x.*((2.*x)./c - 1))./(50.*c.^3.*(1 - (4.*x.^2)./c.^2).^(1./2)) - (77.*t.*(1 - (4.*x.^2)./c.^2).^(1./2))./(100.*c.^2));
theta_u = theta - alpha;
theta_l = -theta - alpha;
cpu = 2.*theta_u / sqrt(M^2-1);
cpl = 2.*theta_l / sqrt(M^2-1);
end

