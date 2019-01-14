set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

% clear all;
clc; close all

select_method = 1;
% select 1 - Newton iteration
% select 2 - Method of lines

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

gamma = 1;
alpha = 2;


% Deposition profile (spatial function for the input)
Pdep1 = (alpha + 1)*gamma*(-(pi/2/L)*sin(pi/2/L*rho)).^alpha.*((pi/2/L)^2*cos(pi/2/L*rho));
Pdep2 = -3/2*cos(pi/2/L*rho);
figure
plot(Pdep1+Pdep2)

% T0 = T0 + Pdep;
T0 = T0 + cos(pi/2/L*rho);
Told = T0;
T = zeros(n+1,m);
T(:,1) = T0;

switch select_method
    case 1
        % Time-varying component for the input
        u1 = exp(-(alpha + 1)*t);
        u2 = exp(-t);
        % Loop for running the Newton algorithm for each time step
        tic
        for j = 2:m % loop thorugh time steps
            T(:,j) = Newton_test(Told,c,drho,dt,Pdep1,Pdep2,u1(j),u2(j),gamma,alpha);
            Told = T(:,j);
        end
        toc
    case 2
        % Time-varying component for the input
        u1 = @(t) exp(-(alpha + 1)*t);
        u2 = @(t) exp(-t);
        % ODE solver for the Method of Lines
        reltol = 1.0e-10; abstol=1.0e-10;
        options = odeset('RelTol',reltol,'AbsTol',abstol);
        tic
        [~,T] = ode23(@ML_test,t,T0,options,drho,n,Pdep1,Pdep2,u1,u2,gamma,alpha);
        toc
        T = T'; % time on collumns
end

% 3D plot
figure
h = surf(t,rho,T)
set(h,'LineStyle','none')
xlabel('t')
ylabel('$\rho$')
zlabel('$T\left(\rho,t\right)$')

Tsol = cos(pi/2/L*rho).*exp(-t) + c;
figure
plot(t,T(1,:)); hold on
plot(t,Tsol(1,:),'r'); hold off

% Err = zeros(length(rho),length(t));

% for i = 1:length(t)
%     for j = 1:length(rho)
%         Err(j,i) = (T(j,i) - Tsol(j,i))^2;
%     end
%     Err(:,i) = Err(:,i) / length(rho);
% end
Err = zeros(1,length(t));
for i = 1:length(t)
   Err(i) = abs(T(1,i) - Tsol(1,i));
end
figure
plot(t,Err);

%%
figure
semilogy(t,[Err_ML; Err_BE]); 
legend('Method of lines','Newton iteration');
title('Error');
xlabel('t [s]');
ylabel('$|T_{sol}-T_{sim}|$');