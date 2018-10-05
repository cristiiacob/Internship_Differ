clc; clear all; close all
L = 50; % Spatial length
t_start = 0.0; t_end = 1;
rho_start = 0; rho_end = L; % 1200 for discrete time
Time = 1; % Time length
a = 0; c = 0; % boundary conditions
n = 999; % spaces in spatial length (grid points = n + 1)
m = 1000; % timesteps

t = linspace(t_start,t_end,m);
rho = linspace(rho_start,rho_end,n+1)'; %TODO check how you define the grid spacing
drho = rho(2) - rho(1);
dt = t(2) - t(1); Fs = 1/dt;
T0 = zeros(n+1,1);
T0(end) = c;
T = zeros(n+1,m);

s = 2/3 * dt / (drho)^3
%%
figure(1)
plot(T0)

% Input
sigma = 5; mu = L/2; K = 10;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
figure(2)
plot(Pdep)

u = sin(2*pi*2*t) + 0.3*sin(2*pi*7*t) + 0.5*sin(2*pi*14*t);

% Initialize Tmperature elements
Told = T0;
T(:,1) = T0;

for j = 2:m % loop thorugh time steps
    T(:,j) = Newton_alg(Told,a,c,drho,dt,Pdep,u(j));
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
S = T(100,:);
Y = fft(S);
P2 = abs(Y/m);
P1 = P2(1:m/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(m/2))/m;
figure(4)
stem(f(1:20),P1(1:20)) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')