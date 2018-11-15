clc; clear all; close all
L = 1000; % length in spatial coordinates
rho_start = 0; rho_end = L;
t_start = 0; t_end = 30; % length in time
a = 0;
C = 2; % constant temperature at the end of the bar
n = 99;  % spaces in length (grid points = n+1)
m = 3000.; % timestamps

t = linspace(t_start,t_end - t_end/m ,m);
rho = linspace(rho_start,rho_end, n+1)';
drho = rho(2) - rho(1);
dt = t(2) - t(1);
T0 = 0*sin(8*pi*rho) + C;
T0(end) = C;

% Pdep = 200*normpdf(0,-round(n/2)+1:round(n/2),50)';
% Ptot = 20; sigma = 8;  MW2keVs = .5;%6.24e21/2.1e19; % Factor to convert to keV
% Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(rho-L/2).^2/sigma.^2);
% sigma = 10; mu = L/2; K = 10;
% sigma = .06; mu = L/2; K = 1/10;
sigma = 100; mu = L/2; K = 10;
Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
plot(Pdep)
u = sin(2*pi*7*t) + sin(2*pi*9*t) + 3; %.* (t<=0.25);
T0 = T0+Pdep;
figure(1)
plot (rho,T0,'-')
title ('Initial temperature distribution')
xlabel ('$\rho$')
ylabel ('T')

options = optimoptions(@fsolve,'Display','iter','MaxIterations',5000,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-5);
Told = T0;
T(:,1) = T0;

for i = 2:2
    [Ti,F,exitflag,output,JAC] = fsolve(@linear_fsolve,T0,options,Told,C,n,drho,dt,Pdep,u(i));
    %     Ti = Ti + Pdep*u(i);
    %     Told = Ti;
    %     T0 = Ti;
    %     T(:,i) = Ti;
end
T = Ti;
% figure(2);
% h = surf(t,rho,T)
% set(h,'LineStyle','none')
%   xlabel('t ')
%   ylabel('\rho')
%   zlabel('T(x,t)')

%%
s = 1/3 * dt / (drho^3);
%    main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
% %    main(1) = -1;
%    main(1) = 2*s*(T(1) - T(2))-1;
%    main(n+1) = 1;
%    upper(1) = 2*s*(T(2) - T(1));
% 
%    for i = 2:n
%       main(i) = 2*s*(T(i-1) - T(i+1))-1;
%       upper(i) = 2*s*(T(i+1) - T(i));
%       under(i-1) = 2*s*(T(i) - T(i-1));
%    end
%    Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
%    Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
%    Under = sparse(2:n+1,1:n,under,n+1,n+1);
%    J = Main + Upper + Under;
%    J = full(J);
%    vals = eig(J)
%    plot(T)


% s = 1000000 / 4 / (drho)^4;
% F(1) = s*(T(2)^3 - 3*T(1)*T(2)^2 + 3*T(1)^2*T(2) - T(1)^3);
% F(n+1) = 0;
% 
% main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
% main(1) = s*(-3*T(2)^2 + 6*T(1)*T(2) - 3*T(1)^2);
% main(n+1) = 0;
% upper(1) = s*(3*T(2)^2 - 6*T(1)*T(2) + 3*T(1)^2);
% 
% for i = 2:n
%     F(i) = s*(T(i+1)^3 - 2*T(i)*T(i+1)^2 - T(i-1)*T(i+1)^2 + 4*T(i-1)*T(i)*T(i+1) - (T(i-1)^2)*T(i+1) - 2*(T(i-1)^2)*T(i) + T(i-1)^3);
%     main(i) = s*(-2*T(i+1)^2 + 4*T(i-1)*T(i+1) - 2*T(i-1)^2);
%     upper(i) = s*(3*T(i+1)^2 - 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) + 4*T(i-1)*T(i) - T(i-1)^2);
%     under(i-1) = s*(-T(i+1)^2 + 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) - 4*T(i-1)*T(i) + 3*T(i-1)^2);
% end
% Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
% Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
% Under = sparse(2:n+1,1:n,under,n+1,n+1);
% J = Main + Upper + Under;
% J = full(J);
% vals = eig(J)
% plot(T)
%
%     s = 1000000  * dt/ 4 / (drho)^4;
%     q = 2/3;
%     n = length(T) - 1;
%     F = zeros(n+1,1);
%     % T = Told;
%     F(1) = s*(T(2)^3 - 3*T(1)*T(2)^2 + 3*T(1)^2*T(2) - T(1)^3) - T(1) + Told(1);
%     F(n+1) = T(n+1) - Told(n+1);
%     
%     main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
%     main(1) = s*(-3*T(2)^2 + 6*T(1)*T(2) - 3*T(1)^2) - 1;
%     main(n+1) = 1;
%     upper(1) = s*(3*T(2)^2 - 6*T(1)*T(2) + 3*T(1)^2);
%     
%     for i = 2:n
%         F(i) = s*(T(i+1)^3 - 2*T(i)*T(i+1)^2 - T(i-1)*T(i+1)^2 + 4*T(i-1)*T(i)*T(i+1) - (T(i-1)^2)*T(i+1) - 2*(T(i-1)^2)*T(i) + T(i-1)^3) - T(i) + Told(i);
%         main(i) = s*(-2*T(i+1)^2 + 4*T(i-1)*T(i+1) - 2*T(i-1)^2) - 1;
%         upper(i) = s*(3*T(i+1)^2 - 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) + 4*T(i-1)*T(i) - T(i-1)^2);
%         under(i-1) = s*(-T(i+1)^2 + 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) - 4*T(i-1)*T(i) + 3*T(i-1)^2);
%     end
%     Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
%     Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
%     Under = sparse(2:n+1,1:n,under,n+1,n+1);
%     J = Main + Upper + Under;
% gamma = 1;
% s = gamma * 1/3 * dt / (drho^3);
% q = 2/3 * dt;
% F = zeros(n+1,1);
% F(1) = s*T(2)  - (s)*T(1);
%     F(n+1) =  T(n+1) - C;
%     
%     main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
%     main(1) = -(s + 1);
%     main(n+1) = 1;
%     upper(1) = s;
%     
%     for i = 2:n
%         F(i) = s*T(i+1) + s*T(i-1) - (2*s)*T(i);
%         main(i) = -(2*s);
%         upper(i) = s;
%         under(i-1) = s;
%     end
%     Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
%     Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
%     Under = sparse(2:n+1,1:n,under,n+1,n+1);
%     J = Main + Upper + Under;

F(1) = s*(T(2)^2 + T(1)^2 - 2*T(1)*T(2));
for i = 2:n
    F(i) = s*(T(i+1)^2 - T(i-1)^2 - 2*T(i)*T(i+1) + 2*T(i-1)*T(i));
end
F(n+1) = T(n+1) - C;

% Evaluate the Jacobian if nargout > 1
    main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
    %    main(1) = -1;
    main(1) = 2*s*(T(1) - T(2));
    main(n+1) = 1;
    upper(1) = 2*s*(T(2) - T(1));
    
    for i = 2:n
        main(i) = 2*s*(T(i-1) - T(i+1));
        upper(i) = 2*s*(T(i+1) - T(i));
        under(i-1) = 2*s*(T(i) - T(i-1));
    end
    Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
    Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
    Under = sparse(2:n+1,1:n,under,n+1,n+1);
    J = Main + Upper + Under;
    
J = full(J);
vals = eig(J)
plot(T)

% ????????????