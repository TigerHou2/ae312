%% Problem 3

close all
clear;clc

disp('Problem 3')
disp(' ')

M = 3;
p = 101e3;
gamma = 1.4;

c = 1; % m
theta  = deg2rad(10);
theta1 = deg2rad(10);
theta2 = deg2rad(-20);

% ------------------------------
% Part a

% calculate top leading edge
beta1 = TBM_get_beta(M,theta1,gamma);
Mn = M * sin(beta1);
M1n = sqrt( (Mn^2 + 2/(gamma-1)) / (2*gamma/(gamma-1)*Mn^2-1) );
M1 = M1n / sin(beta1-theta1);
p1r = TBM_get_pres(M,beta1,gamma);
p1 = p * p1r;
% calculate top trailing edge
mu = asin(1/M1);
fun = @(M2) 0 - theta2 - PM(M2,gamma) + PM(M1,gamma);
M2 = fzero(fun, [1,20]);
p2 = p1 * ( (1+(gamma-1)/2*M1^2) ...
           /(1+(gamma-1)/2*M2^2) ) ^ (gamma/(gamma-1));
% calculate rear end
beta5 = TBM_get_beta(M2,theta1,gamma);
p5r = TBM_get_pres(M2,beta5,gamma);
p5 = p2 * p5r;

disp(['a)  Top/Bottom LE: ' num2str(p1/1e3) ' kPa' newline ...
      '    Top/Bottom TE: ' num2str(p2/1e3) ' kPa' newline ...
      '   Behind Airfoil: ' num2str(p5/1e3) ' kPa' newline ...
      '   Before Airfoil: ' num2str(p /1e3) ' kPa' newline])

    % ------------------------------
    % Part b

    alpha = 0;
    p3 = p1;
    p4 = p2;
    panel = c/2/cos(theta);
    cd = (  sin(theta-alpha)*p1 + sin(-theta-alpha)*p2 ...
           +sin(theta+alpha)*p3 + sin(-theta+alpha)*p4 ) ...
         / (0.5*gamma*p*M^2);
    disp(['b) cd = ' num2str(cd) newline])

% ------------------------------
% Part c

disp(['c) Drag is created due to the pressure differential' newline ...
      '   between the leading the trailing edge of the' newline ...
      '   airfoil. The difference in pressure is caused' newline ...
      '   by the air being accelerated by the expansion' newline ...
      '   fan, which increases the Mach number while' newline ...
      '   decreasing pressure.' newline])
  
    % ------------------------------
    % Part d
    
    % line for oblique shock
    ob = @(x) tan(beta1) * x;
    % line for expansion fan
    xp = @(x) tan(mu+theta) * (x-c/2) + tan(theta)*c/2;
    % intercept
    fun = @(x) ob(x) - xp(x);
    x = fzero(fun,[0,2]);
    y = ob(x);
    
    disp(['d) (x,y) = (' num2str(x) ',' num2str(y) ')' newline])
    
% ------------------------------
% Part e

disp(['e) The shock will bend to the right. Information' newline ...
      '   carried by Mach wave cannnot travel back through' newline ...
      '   the shock because it violates the second law of' newline ...
      '   thermodynamics, and therefore the shock must' newline ...
      '   to the right.'])