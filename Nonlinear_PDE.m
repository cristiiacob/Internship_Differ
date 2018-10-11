clear all; clc; close all
global ncall;

% Parameter initialization
t_start = 0.0; t_end = 1000;
ro_start = 0; ro_end = 50; % 1200 for discrete time
m = 100.; % timestamps
n = 199; % grid spaces ! for Pdep to work, choose odd number
tsim = linspace(t_start,t_end,m);
dt = tsim(2)-tsim(1); Fs = 1/dt;
tspan = [t_start t_end];
ro = linspace(ro_start,ro_end,n+1);
L = ro_end - ro_start; % length of space
dx = L/(n + 1.); % dt = tsim/m;

% Boundary conditions
b = 0; % left boundary condition value
c = 5; % right boundary condition value

% Set input
% Pdep = 100*normpdf(0,-round(n/2)+1:round(n/2),10)'; % deposition profile
sigma = 5; mu = L/2; K = .1;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(ro-mu).^2/sigma.^2)';
% u = @(t) 10*(sin(10*pi*t) + 1)*(t<2.5);
% Ptot = 0.7; sigma = 0.05;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV
% Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(ro-L/2).^2/sigma.^2);
figure(1)
plot(Pdep)
% u = @(t) sin(2*pi*2*t) - sin(2*pi*7*t); % input over time
u = @(t) 1*(1*sin(2*pi*.008*t) + 1) .* (t <= t_end/3);

% Initial profile
T0 = 0*sin(2*pi/100*ro) + c;
T0(n+1) = c; % constant boundary condition (impose)

% ODE call
reltol = 1.0e-10; abstol=1.0e-10;
options = odeset('RelTol',reltol,'AbsTol',abstol);
% tsim is used in order to force fixed time step
% [t,T] = ode45(@nonlin_eq,tsim,T0,options,u,Pdep,b,dx,n);
[t,T] = ode15s(@lin_eq,tsim,T0,options,u,Pdep,b,dx,n);
% [t,T] = ode45(@nonlin_system,tsim,T0,options,u,Pdep,b,dx,n);
T = T'; % time on collumns
figure(2);
h = surf(t,ro,T)
set(h,'LineStyle','none')
xlabel('t ')
ylabel('\rho')
zlabel('T(x,t)')

%%               
S = T(100,:);
Y = fft(S);
P2 = abs(Y/m);
P1 = P2(1:m/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(m/2))/m;
figure(3)
stem(f(1:10),P1(1:10)) 
set(gca,'YScale','log');
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%%
% clear all; clc; close all
% 
% % Parameter initialization
% t_start = 0.0; t_end = 1;
% ro_start = 0; ro_end = 100; % 1200 for discrete time
% m = 100.; % timestamps
% n = 199; % grid spaces ! for Pdep to work, choose odd number
% tsim = linspace(t_start,t_end,m);
% dt = tsim(2)-tsim(1); Fs = 1/dt;
% tspan = [t_start t_end];
% ro = linspace(ro_start,ro_end,n+1);
% L = ro_end - ro_start; % length of space
% dx = L/(n + 1.); % dt = tsim/m;
% 
% % Boundary conditions
% b = 0; % left boundary condition value
% c = 1; % right boundary condition value
% 
% % Set input
% % Pdep = 100*normpdf(0,-round(n/2)+1:round(n/2),10)'; % deposition profile
% sigma = 10; mu = L/2; K = 101;
% Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(ro-mu).^2/sigma.^2)';
% % u = @(t) 10*(sin(10*pi*t) + 1)*(t<2.5);
% % Ptot = 0.7; sigma = 0.05;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV
% % Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(ro-L/2).^2/sigma.^2);
% figure(1)
% plot(Pdep)
% % u = @(t) sin(2*pi*2*t) - sin(2*pi*7*t); % input over time
% u = @(t) (1*sin(2*pi*5*t)) %.* (t <= t_end/3);
% 
% % Initial profile
% T0 = 0*sin(2*pi/100*ro) + c + Pdep'*u(1);
% T0(n+1) = c; % constant boundary condition (impose)
% 
% % ODE call
% reltol = 1.0e-10; abstol=1.0e-10;
% options = odeset('RelTol',reltol,'AbsTol',abstol);
% % tsim is used in order to force fixed time step
% % [t,T] = ode45(@nonlin_eq,tsim,T0,options,u,Pdep,b,dx,n);
% % [t,T] = ode15s(@lin_eq,tsim,T0,options,u,Pdep,b,dx,n);
% [t,T] = ode45(@nonlin_system,tsim,T0,options,u,Pdep,b,dx,n);
% T = T'; % time on collumns
% figure(2);
% h = surf(t,ro,T)
% set(h,'LineStyle','none')
% xlabel('t ')
% ylabel('\rho')
% zlabel('T(x,t)')