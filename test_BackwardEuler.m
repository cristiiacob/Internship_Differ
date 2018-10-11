clc; clear all; close all
L = .5; % length in spatial coordinates
ro_start = 0; ro_end = L; 
Time = .5; % length in time
a = 0;  
C = 5; % constant temperature at the end of the bar
n = 499;  % grid points in space
m = 1000.; % timestamps

dx = L/(n + 1.); dt = Time/m; 

t = 0.;
P = 200*normpdf(0,-round(n/2):round(n/2),50)';
% pvec = 150*normpdf(0,-m/2+1:m/2,75);
time_v = 0:m-1;
pvec = 0*.1*sin(.01*pi*time_v);
for j = 1:n+2
   pvect(j,:) = pvec;  
end
gamma = 1; % Diffusivity constant

q = gamma * dt / (dx^2) * 2/3; 

Adiag = (2*q + 1)*ones(n+2,1); Adiag(end) = 1;
Aover = -q*ones(n+2,1); Aover(2) = -2.*q;
Aunder = -q*ones(n+2,1); Aunder(n+1:end) = 0;

A = spdiags([Aunder,Adiag,Aover],[-1 0 1],n+2,n+2);

b = zeros(n+2,1);
b(1) = -2.*q*a*dx;
% b(end) = q*C;

ro = linspace(ro_start,ro_end, n+2)';

% T_tk = 1*sin(2*pi*ro) + 5;
T_tk = 100*normpdf(0,-round(n/2):round(n/2),50)' + 5;
T_tk(end) = C;

figure(1)
plot (ro,T_tk,'-')
title ('Initial condition for Temperature distribution')
xlabel ('ro')
ylabel ('T')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=zeros(n+2,m);
tvec=zeros(m,1);

T(:,1) = T_tk;
tvec(1) = t;

for k = 1:m
    t = t+dt;
    % if boundary conditions vary with time you need to update b here
    % with implicit method we solve a matrix equation at each step:
    c = T_tk + b;
    T_tk_1 = A\c;

        T_tk_1 = T_tk_1 + P.*pvect(:,k);

    % this is very time consuming later we will discuss faster ways to solve this problem using iterative methods

    T(:,k) = T_tk_1;
 
    % for next time step:
    T_tk = T_tk_1;
    tvec(k) = t;

end

figure(2)
mesh (tvec,ro,T)
title ('Temperature distribution over time')
xlabel ('t')
ylabel ('rho')
zlabel ('T')