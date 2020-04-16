%% Problem 3

close all
clear;clc

syms x real

t = 0.03;
c = 1;
Mi = 2.5;
p = 100;
gamma = 1.4;

h = t - t.*(x./c).^2;
theta = diff(h,x);
theta0 = atan(double(subs(theta,x,-1))); % theta at the leading edge

pres_points = 25;
syms m1 real
mU = sym('mU', [pres_points,1]) ;
mL = sym('mL', [pres_points,1]) ;

alpha_vect = deg2rad(linspace(0,12,25));
cl_vect = zeros(length(alpha_vect),1);
cd_vect = zeros(length(alpha_vect),1);

% alpha_vect = deg2rad([0.1]);

for i = 1:length(alpha_vect)
    
    alpha = alpha_vect(i);

% top surface LE oblique shock / expansion fan
if alpha <= theta0 % oblique shock
    t1 = theta0 - alpha;
    beta1 = TBM_get_beta(Mi,t1,gamma);
    p1 = p * TBM_get_pres(Mi,beta1,gamma);
    Min = Mi * sin(beta1);
    M1n = sqrt( (Min^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*Min^2-1) );
    M1 = M1n / sin(beta1-t1);
else % alpha > theta0, expansion fan
    t1 = theta0 - alpha;
    eqn = -t1 == PM(m1,gamma) - PM(Mi,gamma);
    M1 = abs(double(vpasolve(eqn,m1)));
    p1 = p * ( (1+(gamma-1)/2*Mi^2) ...
              /(1+(gamma-1)/2*M1^2) ) ^ (gamma/(gamma-1));
end

% bottom surface LE oblique shock
    t3 = theta0 + alpha;
    beta3 = TBM_get_beta(Mi,t3,gamma);
    p3 = p * TBM_get_pres(Mi,beta3,gamma);
    Min = Mi * sin(beta3);
    M3n = sqrt( (Min^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*Min^2-1) );
    M3 = M3n / sin(beta3-t3);

% sample a whole bunch of points for pressure distribution
x = linspace(-1,1,pres_points)';
theta_vect = double(subs(theta,x));

thetas = theta_vect-theta_vect(1);
eqnU = -thetas == PM(mU,gamma) - PM(M1,gamma);
eqnL = -thetas == PM(mL,gamma) - PM(M3,gamma);

MU = struct2cell(vpasolve(eqnU,mU));
ML = struct2cell(vpasolve(eqnL,mL));

m_uu = zeros(length(MU),1);
m_ll = zeros(length(ML),1);

for k = 1:length(MU)
    m_uu(k) = abs(double(MU{k}));
    m_ll(k) = abs(double(ML{k}));
end

pu = p1 .* ( (1+(gamma-1)./2.*  M1.^2) ...
           ./(1+(gamma-1)./2.*m_uu.^2) ) .^ (gamma./(gamma-1));
pl = p3 .* ( (1+(gamma-1)./2.*  M3.^2) ...
           ./(1+(gamma-1)./2.*m_ll.^2) ) .^ (gamma./(gamma-1));

cl = sum( pl.*cos(theta_vect+alpha)...
         -pu.*cos(theta_vect-alpha) )...
    /length(pu)/(0.5*gamma*p*Mi^2);
cd = sum( pl.*sin(theta_vect+alpha)...
         +pu.*sin(theta_vect-alpha) )...
    /length(pu)/(0.5*gamma*p*Mi^2);

cl_vect(i) = cl;
cd_vect(i) = cd;

end

figure(3)
plot(rad2deg(alpha_vect),cl_vect./cd_vect)
title('$C_L/C_D$ vs. $\alpha$')
xlabel('$\alpha$')
ylabel('$C_L/C_D$')
grid(gca,'minor')
grid on
latexify