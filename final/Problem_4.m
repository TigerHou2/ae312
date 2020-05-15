%% Problem 4

close all
clear;clc

disp('Problem 4')
disp(' ')

gamma = 1.4;
R  = 287;

At = 37.93e-4; % m^2
Ae = 100e-4;   % m^2

p0 = 800e3; % Pa

% ------------------------------
% Part a

gamm = (gamma-1)/2;
gamp = (gamma+1)/2;

fun = @(M) 1/M^2 * (1/gamp * (1+gamm * M^2))^(gamp/gamm) - (Ae/At)^2;
Me3 = fzero(fun, [1e-16,1]);
Me6 = fzero(fun, [1,20]);
Pb3 = (1+gamm*Me3^2)^(-gamma/2/gamm);
Pb6 = (1+gamm*Me6^2)^(-gamma/2/gamm);
Pbmin = Pb3*p0;
disp(['a) Pb_min = ' num2str(Pbmin/1000) ' kPa' newline])

    % ------------------------------
    % Part b
    Pe = Pb6 * p0;
    Pb = Pe;
    Me = Me6;
    disp(['b) Pe = ' num2str(Pe/1e3) ' kPa' newline ...
          '   Pb = ' num2str(Pb/1e3) ' kPa' newline ...
          '   Me = ' num2str(Me) newline])

% ------------------------------
% Part c

% we know from part a that the flow is subsonic for 780 kPa
Pb = 780e3;
Pe = Pb;
Me = sqrt(1/gamm * ((p0/Pe)^((gamma-1)/gamma)-1));
% mass flow rate is the same at throat and exit
fun = @(Mt) At*p0*sqrt(gamma/R)*Mt/((1+gamm*Mt^2)^(gamp/gamm/2))...
          - Ae*p0*sqrt(gamma/R)*Me/((1+gamm*Me^2)^(gamp/gamm/2));
Mt = fzero(fun,[1e-16,1]);
disp(['c) Me = ' num2str(Me) newline ...
      '   Mt = ' num2str(Mt) newline])
  
    % ------------------------------
    % Part d

    Me = Me6;
    Te = 283;
    % no shock at design point, can use isentropic chain
    T0 = Te * (1+gamm*Me^2);
    mdot = Ae*p0/sqrt(T0)*sqrt(gamma/R)*Me/((1+gamm*Me^2)^(gamp/gamm/2));
    disp(['d) T0   = ' num2str(T0) ' K' newline ...
          '   mdot = ' num2str(mdot) ' kg/s' newline])
  
% ------------------------------
% Part e
p0vect = [1200e3 800e3  600e3  300e3  100e3  10e3];
mdvect = [7.2907 4.8605 3.6454 1.8227 0.6076 0.0608];
plot(p0vect,mdvect)

    % ------------------------------
    % Part f
    Me5 = sqrt((Me6^2 + 1/gamma) / (gamma/gamm * Me6^2 - 1));
    Pb5 = Pb6 * (1+gamma/gamp*(Me6^2-1));
    Pb = 500e3;
    if Pb/p0 >= Pb5 && Pb/p0 < Pb3
        disp('e) Yes, there is a shock for Pb = 500 kPa')
    else
        disp('e) No, there is no shock for Pb = 500 kPa')
    end
    Pb = 200e3;
    if Pb/p0 >= Pb5 && Pb/p0 < Pb3
        disp('   Yes, there is a shock for Pb = 200 kPa')
    else
        disp('   No, there is no shock for Pb = 200 kPa')
    end
    disp(' ')
    
