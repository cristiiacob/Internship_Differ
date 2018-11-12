clear all; clc; close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex'); 
global ncall;

% Parameter initialization
t_start = 0.0; t_end = 10;
ro_start = 0; ro_end = .5; % 1200 for discrete time
m = 1000.; % timestamps
n = 499; % grid spaces ! for Pdep to work, choose odd number
tsim = linspace(t_start,t_end - t_end/m,m);
dt = tsim(2)-tsim(1); Fs = 1/dt;
tspan = [t_start t_end];
ro = linspace(ro_start,ro_end,n+1);
L = ro_end - ro_start; % length of space
dx = L/(n + 1.); % dt = tsim/m;

% Boundary conditions
b = 0; % left boundary condition value
c = 2; % right boundary condition value

% Set input
% Pdep = 100*normpdf(0,-round(n/2)+1:round(n/2),10)'; % deposition profile
sigma = .06; mu = L/2; K = 1;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(ro-mu).^2/sigma.^2)';
% Pdep = 200*normpdf(0,-round(n/2):round(n/2),50)';
% u = @(t) 10*(sin(10*pi*t) + 1)*(t<2.5);
% Ptot = 0.7; sigma = 0.05;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV
% Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(ro-L/2).^2/sigma.^2);
figure(1)
plot(Pdep)
% u = @(t) sin(2*pi*2*t) - sin(2*pi*7*t); % input over time
u = @(t) 1*(sin(2*pi*7*t) + sin(2*pi*9*t) + 3);

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
ylabel('$\rho$')
zlabel('$T\left(\rho,t\right)$')

%%               
 m = 200;
S = T((n+1)/2,801:end);
Y = fft(S)/length(S);
P2 = abs(Y);
P1 = P2(1:m/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(m/2))/m;
figure(4)

stem(f,P1) ;
set(gca,'YScale','log'); 
title('Single-Sided Amplitude Spectrum of T(t)')
xlabel('f(Hz)')
ylabel('$|T(f)|$')

figure(5)
plot(T((n+1)/2,:))

%% 2nd order
% clear all; clc; close all
% 
% % Parameter initialization
% t_start = 0.0; t_end = 10;
% ro_start = 0; ro_end = 1000; % 1200 for discrete time
% m = 1000.; % timestamps
% n = 199; % grid spaces ! for Pdep to work, choose odd number
% tsim = linspace(t_start,t_end - t_end/m,m);
% dt = tsim(2)-tsim(1); Fs = 1/dt;
% tspan = [t_start t_end];
% ro = linspace(ro_start,ro_end,n+1);
% L = ro_end - ro_start; % length of space
% dx = L/(n + 1.); % dt = tsim/m;
% % Boundary conditions
% b = 0; % left boundary condition value
% c = 2; % right boundary condition value
% 
% % Set input
% % Pdep = 100*normpdf(0,-round(n/2)+1:round(n/2),10)'; % deposition profile
% sigma = 100; mu = L/2; K = 10000;
% Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(ro-mu).^2/sigma.^2)';
% % u = @(t) 10*(sin(10*pi*t) + 1)*(t<2.5);
% % Ptot = 0.7; sigma = 0.05;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV
% % Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(ro-L/2).^2/sigma.^2);
% figure(1)
% plot(Pdep)
% 
% % u = @(t) sin(2*pi*2*t) - sin(2*pi*7*t); % input over time
% u = @(t) 1*(sin(2*pi*7*t) + sin(2*pi*9*t) + 3);
% 
% % Initial profile
% T0 = 0*sin(2*pi/100*ro) + c;
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
% ylabel('$\rho$')
% zlabel('$T\left(\rho,t\right)$')

%% 3rd order
clear all; clc; close all

% Parameter initialization
t_start = 0.0; t_end = 50;
ro_start = 0; ro_end = 1000; % 1200 for discrete time
m = 5000.; % timestamps
n = 99; % grid spaces ! for Pdep to work, choose odd number
tsim = linspace(t_start,t_end - t_end/m,m);
dt = tsim(2)-tsim(1); Fs = 1/dt;
tspan = [t_start t_end];
ro = linspace(ro_start,ro_end,n+1);
L = ro_end - ro_start; % length of space
dx = L/(n + 1.); % dt = tsim/m;
% Boundary conditions
b = 0; % left boundary condition value
c = 2; % right boundary condition value

% Set input
% Pdep = 100*normpdf(0,-round(n/2)+1:round(n/2),10)'; % deposition profile
sigma = 125; mu = L/2; K = 100000;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(ro-mu).^2/sigma.^2)'; Pdep(1) = 0; Pdep(end) = 0;
% u = @(t) 10*(sin(10*pi*t) + 1)*(t<2.5);
% Ptot = 0.7; sigma = 0.05;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV
% Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(ro-L/2).^2/sigma.^2);
figure(1)
plot(Pdep)

% u = @(t) sin(2*pi*2*t) - sin(2*pi*7*t); % input over time
u = @(t) 1*(sin(2*pi*7*t) + sin(2*pi*9*t) + 3);

% Initial profile
T0 = zeros(n+1,1) + c;
% T0 = T0 + Pdep;
T0(1) = c; T0(n+1) = c; % constant boundary condition (impose)

% ODE call
reltol = 1.0e-10; abstol=1.0e-10;
options = odeset('RelTol',reltol,'AbsTol',abstol);
% tsim is used in order to force fixed time step
% [t,T] = ode45(@nonlin_eq,tsim,T0,options,u,Pdep,b,dx,n);
% [t,T] = ode15s(@lin_eq,tsim,T0,options,u,Pdep,b,dx,n);
[t,T] = ode45(@nonlin_sys_3rd,tsim,T0,options,u,Pdep,b,dx,n);
T = T'; % time on collumns
figure(2);
h = surf(t,ro,T)
set(h,'LineStyle','none')
xlabel('t ')
ylabel('$\rho$')
zlabel('$T\left(\rho,t\right)$')
figure(3)
plot(T((n+1)/2,:));

%%               
 m = 600;
S = T((n+1)/2-20,4401:end);
Y = fft(S)/length(S);
P2 = abs(Y);
P1 = P2(1:m/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(m/2))/m;
figure(4)

stem(f,P1) ;
set(gca,'YScale','log'); 
title('Single-Sided Amplitude Spectrum of T(t)')
xlabel('f(Hz)')
ylabel('$|T(f)|$')

figure(5)
plot(T((n+1)/2-20,:))