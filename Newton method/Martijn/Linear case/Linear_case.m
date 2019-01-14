set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

clear all; clc; close all

L = 1; % Spatial length
t_start = 0.0; t_end = 10; % time limits
rho_start = 0; rho_end = L; % spatial limits
c = 2; % constant boundary condition
n = 99; % spaces in spatial length (grid points = n + 1)
m = 10000; % time points

t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end,n+1)';
drho = rho(2) - rho(1);
dt = t(2) - t(1); Fs = 1/dt;
T0 = zeros(n+1,1) + c; % initial temperature distribtution
Told = T0;
T = zeros(n+1,m);
T(:,1) = T0;
% Deposition profile (spatial function for the input)
sigma = .125; mu = L/2; K = 100000;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2); Pdep(1) = 0; Pdep(end) = 0;
% figure
% plot(Pdep)

%% Newton method
        % Time-varying component for the input
        u = sin(2*pi*11*t) + 1*sin(2*pi*14*t) + 3;
        
        % Loop for running the Newton algorithm for each time step
        tic
        for j = 2:m % loop thorugh time steps
            T(:,j) = Newton_linear(Told,c,drho,dt,Pdep,u(j));
            Told = T(:,j);
        end
        toc

%% Plots

% 3D plot
figure
h = surf(t,rho,T)
set(h,'LineStyle','none')
xlabel('t')
ylabel('$\rho$')
zlabel('$T\left(\rho,t\right)$')
title('Temperature distribution evolution')

% figure
% 
% subplot(121)
% % time plot of the temperature evolution of one spatial point
% plot(t,T((n+1)/2-20,:))
% xlabel('t')
% ylabel('$T\left(30,t\right)$')
% 
% subplot(122)
% % spatial temperature distribution at a particular time instant
% plot(rho,T(:,600))
% xlabel('$\rho$')
% ylabel('$T\left(\rho,t\right)$')