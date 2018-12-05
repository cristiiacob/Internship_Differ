clc; clear all; close all
L = 1; % Spatial length
t_start = 0.0; t_end = .1;
rho_start = 0; rho_end = L; % 1200 for discrete time
Time = 1; % Time length
a = 0; c = 1; % boundary conditions
n = 99; % spaces in spatial length (grid points = n + 1)
m = 1000; % time points
t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end,n+1)';
r = 20;
C = eye(r+1);
gamma =  1;
C1 = -gamma * (2*pi/L)^2;
diagonal = 0:r;
diagonal = diagonal .* diagonal;
A = C1*diag(diagonal);

phi_k = @(rho,n) (n>=1)*(sqrt(2/L)/c*cos(2*n*pi*rho/L)) + (n < 1)*1/c/sqrt(L);
sigma = .075; mu = L/2; K = 1;
Pdep = @(rho) K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2); 
T0_init = @(rho) c + 1*K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
% T0_init = @(rho) cos(2*pi*rho*3/L/4);
% figure(1)
% fplot(Pdep,[rho_start rho_end]);
T0 = T0_init(rho);
figure(1);
plot(rho,T0);
T = zeros(n+1,m);
T(:,1) = T0; Ti = zeros(n+1,m);
a0 = zeros(r+1,1);
B = zeros(r+1,1);

figure(3)
for i = 0:r
    B(i+1) = integral(@(rho) Pdep(rho).*phi_k(rho,i),0,L);
    a0(i+1) = integral(@(rho) T0_init(rho).*phi_k(rho,i),0,L);
    fplot(@(rho) phi_k(rho,i),[rho_start rho_end]); hold on
end
hold off;


Test = zeros(n+1,1); Test_j = Test;
for i = 0:r
    Test_j = a0(i+1)*phi_k(rho,i);
    Test = Test + Test_j;
    Test_j = zeros(n+1,1);
end
figure(4)
plot(rho,Test)
%%
sys = ss(A,B,C,[]);
u = 0*(sin(2*pi*7*t) + sin(2*pi*9*t) + 3);
a = lsim(sys,u,t,a0); a = a';

for i = 0:r
    Ti(:,2:end) = a(i+1,2:end).*phi_k(rho,i);
    T = T + Ti;
    Ti = zeros(n+1,m);
end

figure(4)
h = surf(t,rho,T)
set(h,'LineStyle','none')
xlabel('t ')
ylabel('$\rho$')
zlabel('T(x,t)')
