clc; clear all; close all
L = .5; % length in spatial coordinates
rho_start = 0; rho_end = L;
Time = 1; % length in time
a = 0;
C = 5; % constant temperature at the end of the bar
n = 199;  % spaces in length (grid points = n+1)
m = 100.; % timestamps

t = linspace(0,Time,m);
rho = linspace(rho_start,rho_end, n+1)';
drho = rho(2) - rho(1);
dt = t(2) - t(1);
T0 = 0*sin(8*pi*rho) + 5;
T0(end) = C;

% Pdep = 200*normpdf(0,-round(n/2)+1:round(n/2),50)';
Ptot = 0.7; sigma = 0.05;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV   
Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(rho-L/2).^2/sigma.^2);
plot(Pdep)
u = sin(2*pi*2*t) + 0.3*sin(2*pi*7*t);

figure(1)
plot (rho,T0,'-')
title ('Initial temperature distribution')
xlabel ('\rho')
ylabel ('T')

options = optimoptions(@fsolve,'Display','iter','MaxIterations',300,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-3);
Told = T0;
T(:,1) = T0;

for i = 2:length(t)
    [Ti,F,exitflag,output,JAC] = fsolve(@nl_discretized,T0,options,Told,C,n,drho,dt);
    Ti = Ti + Pdep*u(i);
    Told = Ti;
    T0 = Ti;
    T(:,i) = Ti;
end

figure(2);
h = surf(t,rho,T)
set(h,'LineStyle','none')
  xlabel('t ')
  ylabel('\rho')
  zlabel('T(x,t)')