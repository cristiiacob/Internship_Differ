clc; clear all; close all
% the nonlinear PDE with alpha = 1
L = 1000; % Spatial length
t_start = 0.0; t_end = 10;
rho_start = 0; rho_end = L; % 1200 for discrete time
Time = 1; % Time length
a = 0; c = 2; % boundary conditions
n = 199; % spaces in spatial length (grid points = n + 1)
m = 1000; % time points 

% t = linspace(t_start,t_end ,m);
t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end,n+1)'; %TODO check how you define the grid spacing
drho = rho(2) - rho(1);
dt = t(2) - t(1); Fs = 1/dt;
T0 = zeros(n+1,1) + c;
T0(end) = c;
% T0 = bvp_test';
% T0 = T0(:,1);
T = zeros(n+1,m);

s = 2/3 * dt / (drho)^3
%%
% figure(1)
% plot(T0)

% Input
sigma = 100; mu = L/2; K = 10000;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
figure(1)
plot(Pdep)

u = 0*(sin(2*pi*7*t) + sin(2*pi*9*t) + 3);

% Initialize Tmperature elements
% T0 = T0 + Pdep*u(1) ;
figure(2)
T0 = T0 + Pdep;
plot(T0)
Told = T0;
T(:,1) = T0;

for j = 2:m % loop thorugh time steps
    T(:,j) = Crank_Nicholson(Told,a,c,drho,dt,Pdep,u(j));
    T(:,j) = T(:,j) ; % Should introduce in Newton?
    Told = T(:,j);
end

figure(3);
h = surf(t,rho,T)
set(h,'LineStyle','none')
xlabel('t ')
ylabel('\rho')
zlabel('T(x,t)')

%%               
S = T((n+1)/2,:);
Y = fft(S)/length(S);
P2 = abs(Y);
P1 = P2(1:m/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(m/2))/m;
figure(4)

stem(f,P1) ;
set(gca,'YScale','log'); 
title('Single-Sided Amplitude Spectrum of T(t)')
xlabel('f (Hz)')
ylabel('|T(f)|')

figure(5)
plot(T((n+1)/2,:))
%%
% figure(5)
% t = linspace(0,10,10)
% T = 1;
% Fs = 1/T;
% x = sin(2*pi*t/T) + 1;
% X = fft(x) / length(x);
% X = abs(X);
% stem(X)

