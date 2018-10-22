clc; clear all; close all
L = .5; % length in spatial coordinates
rho_start = 0; rho_end = L;
t_start = 0.0; t_end = 10;
a = 0;  
C = 2; % constant temperature at the end of the bar
n = 498;  % grid points in space
m = 1000.; % timestamps

t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end,n+2)'; %TODO check how you define the grid spacing
drho = rho(2) - rho(1);
dt = t(2) - t(1); Fs = 1/dt;

P = 200*normpdf(0,-round(n/2):round(n/2)+1,50)';
pvec = sin(2*pi*7*t) + sin(2*pi*9*t) + 3;
for j = 1:n+2
   pvect(j,:) = pvec;  
end
gamma = 1; % Diffusivity constant

q = gamma * dt / (drho^2) * 2/3; 

Adiag = (2*q + 1)*ones(n+2,1); Adiag(end) = 1;
Aover = -q*ones(n+2,1); Aover(2) = -2.*q;
Aunder = -q*ones(n+2,1); Aunder(n+1:end) = 0;

A = spdiags([Aunder,Adiag,Aover],[-1 0 1],n+2,n+2);

b = zeros(n+2,1);
b(1) = -2.*q*a*drho;
% b(end) = q*C;

rho = linspace(rho_start,rho_end, n+2)';

T_tk = 0*100*normpdf(0,-round(n/2):round(n/2)+1,50)' + C;
T_tk(end) = C;

figure(1)
plot (rho,T_tk,'-')
title ('Initial condition for Temperature distribution')
xlabel ('ro')
ylabel ('T')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=zeros(n+2,m);
T(:,1) = T_tk;

for k = 1:m
    c = T_tk + b;
    T_tk_1 = A\c;
    T_tk_1 = T_tk_1 + P.*pvect(:,k);
    T(:,k) = T_tk_1;
    T_tk = T_tk_1;
end

figure(2)
mesh (t,rho,T)
title ('Temperature distribution over time')
xlabel ('t')
ylabel ('rho')
zlabel ('T')

figure(3)
plot(T((n+2)/2,:))

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
ylabel('|T(f)|')
