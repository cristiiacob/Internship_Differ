clc; clear all; close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex'); 
L = .5; % length in spatial coordinates
% L = 1;
rho_start = 0; rho_end = L;
t_start = 0.0; t_end = 1;
a = 0;  
C = 2; % constant temperature at the end of the bar
n = 498;  % grid points in space
m = 1000.; % timestamps

t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end,n+2)'; %TODO check how you define the grid spacing
drho = rho(2) - rho(1);
dt = t(2) - t(1); Fs = 1/dt;

% P = 200*normpdf(0,-round(n/2):round(n/2)+1,50)';
sigma = .06; mu = L/2; K = 1;
P = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2); P(1) = 0; P(end) = 0;
% Pdep = ones(1,n+1); Pdep(end) = 0;
figure(1)
plot(P)
pvec = sin(2*pi*7*t) + sin(2*pi*9*t) + 3;
for j = 1:n+2
   pvect(j,:) = 0*pvec;  
end
gamma = 1/2; % Diffusivity constant
% gamma = 1;

q = gamma * dt / (drho^2) * 2/3;

Adiag = (2*q + 1)*ones(n+2,1); Adiag(end) = 1;
Aover = -q*ones(n+2,1); Aover(2) = -2.*q;
Aunder = -q*ones(n+2,1); Aunder(n+1:end) = 0;

A = spdiags([Aunder,Adiag,Aover],[-1 0 1],n+2,n+2);

b = zeros(n+2,1);
b(1) = -1.*q*a*drho;
% b(end) = q*C;

rho = linspace(rho_start,rho_end, n+2)';

T_tk = 0*100*normpdf(0,-round(n/2):round(n/2)+1,50)' + C + P;
T_tk(end) = C;

figure(1)
plot (rho,T_tk,'-')
title ('Initial condition for Temperature distribution')
xlabel ('$\rho$')
ylabel ('T')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=zeros(n+2,m);
T(:,1) = T_tk;

for k = 2:m
    c = T_tk + b + dt*2/3*P.*pvect(:,k-1);
    T_tk_1 = A\c;
    T_tk_1 = T_tk_1 ;
    T(:,k) = T_tk_1;
    T_tk = T_tk_1;
end

figure(2)
mesh (t,rho,T)
title ('Temperature distribution over time')
xlabel ('t')
ylabel ('$\rho$')
zlabel ('$T\left(\rho,t\right)$')

figure(3)
plot(T((n+2)/2,:))
T(1,80)
 m = 200;
S = T((n+2)/2,801:end);
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
ylabel('$|T\left(f\right)|$')

figure(5)
plot(T(1,:))