clc; clear all; close all
L = .5; % Spatial length
rho_start = 0; rho_end = L;
Time = 1; % Time length
a = 0;
c = 5;
n = 199; % spaces in spatial length (grid points = n + 1)
m = 100; % timesteps

t = linspace(0,Time,m);
rho = linspace(rho_start,rho_end,n+1)';
drho = rho(2) - rho(1);
dt = t(2) - t(1);
T0 = 0*sin(8*pi*rho) + 5;
T0(end) = c;
T = zeros(n+1,m);

% Input
sigma = 0.05; mu = L/2; K = .1;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
plot(Pdep)
u = sin(2*pi*2*t) + 0.3*sin(2*pi*7*t);

% Initialize Tmperature elements
Told = T0;
T(:,1) = T0;

for j = 2:m % loop thorugh time steps
    T_temp = Newton_alg(Told,a,b,drho,dt);
    T(:,j) = T_temp;
    Told = T_temp;
end

figure(2);
h = surf(t,rho,T)
set(h,'LineStyle','none')
xlabel('t ')
ylabel('\rho')
zlabel('T(x,t)')