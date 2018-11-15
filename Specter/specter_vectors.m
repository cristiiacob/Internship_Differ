clc; clear all; close all
L = .5; % Spatial length
t_start = 0.0; t_end = 1;
rho_start = 0; rho_end = L; % 1200 for discrete time
Time = 1; % Time length
a = 0; c = 2; % boundary conditions
n = 499; % spaces in spatial length (grid points = n + 1)
m = 1000; % time points
t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end,n+1)';
r = 5;
C = eye(r+1);
gamma =  1e-2;
C1 = -gamma * (2*pi/L)^2;
diagonal = 1:(r+1);
diagonal = diagonal .* diagonal;
A = C1*diag(diagonal);
phi_k = zeros(r+1,length(rho));

for i = 1:r+1
    phi_k(i,:) = (i>1)*(sqrt(2/L)/c*cos(2*i*pi*rho/L)) + (i == 1)*1/sqrt(L);
end
sigma = .05; mu = L/2; K = 1;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
T0_init = c + 0 * K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
% T0_init = @(rho) cos(2*pi*rho*3/L/4);
figure(1)
plot(rho,Pdep);
T0 = T0_init;
figure(2);
plot(rho,T0);
T = zeros(n+1,m);
T(:,1) = T0; Ti = zeros(n+1,m);
a0 = zeros(r+1,1);
B = zeros(r+1,1);

figure(3)
for i = 1:r+1
    B(i) = trapz(rho,Pdep'.*phi_k(i,:));
    a0(i) = trapz(rho,T0_init'.*phi_k(i,:));
    plot(rho,phi_k(i,:)); hold on
end
hold off;


Test = zeros(n+1,1); Test_j = Test;
for i = 1:r
    for j = 1:length(rho)
        Test_j(j) = a0(i)*phi_k(i+1,j);
    end
    Test = Test + Test_j;
    Test_j = zeros(n+1,1);
end
figure(4)
plot(Test)
%%
sys = ss(A,B,C,[]);
u = 0*(sin(2*pi*7*t) + sin(2*pi*9*t) + 3);
a = lsim(sys,u,t,a0); a = a';


for i = 1:r+1
    for j = 1:length(rho)
        Ti(j,2:end) = a(i,2:end)*phi_k(i,j);
    end
    T = T + Ti;
    Ti = zeros(n+1,m);
end

figure(4)
h = surf(t,rho,T)
set(h,'LineStyle','none')
xlabel('t ')
ylabel('\rho')
zlabel('T(x,t)')
