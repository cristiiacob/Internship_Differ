clc; clear all; close all
length = 1;
tmax = 2;
jmax = 2500;
dt = tmax / jmax;
nmax = 250;
dx = length / nmax;
q1 = .05;alpha = q1 * dt / (dx)^2;

grid = linspace(0,length,nmax+1);
time = linspace(0,tmax,jmax+1);
T = zeros(nmax+1,jmax+1);

f1 = 5; f2 = 3;
% Initial conditiom
for n = 1:nmax+1
   T(n,1) = .5*sin(2*pi*f1*grid(n))+sin(2*pi*f2*grid(n))+5;
end
% Boundary condition
T(1,:) = 0;
T(nmax+1,:) = 0;
a(1:nmax-2)=-alpha;
b(1:nmax-1)=1+2*alpha;
c(1:nmax-2)=-alpha;
M=inv(diag(b,0)+diag(a,-1)+diag(c,1));

for k=2:jmax % Time Loop
    t_prev=T(2:nmax,k-1);
    T(2:nmax,k)=M*t_prev;
end

figure(1)
plot(grid,T(:,1),'-',grid,T(:,100),'-',grid,T(:,300),'-',grid,T(:,600),'-')
title('Temperature within the fully implicit method')
xlabel('X')
ylabel('T')
figure(2)
mesh(grid,time,T')
title('Temperature within the fully implicit method')
xlabel('grid')
ylabel('time')