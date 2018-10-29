gradT = -(0:.0001:.02);
chi = 1/2;
alpha = 1/2;

q = (-chi * gradT) .^ (alpha+1);

plot(-gradT,q);
%%
syms gradT_sym
q_sym = (-chi * gradT_sym) ^ (alpha+1);

pol_taylor = taylor(q_sym,gradT_sym,'ExpansionPoint',-0.01,'Order',4)
pol_coeffs = collect(pol_taylor,gradT_sym)
q_fitted = subs(pol_coeffs,gradT_sym,gradT);
plot(-gradT,[q; q_fitted]);



%% test initial temperature profile - not finding solutions
% syms gradT(rho) P alpha C
% eqn = diff(gradT,rho,2) == -P / (alpha+1) * 1/(gradT ^ alpha);
% DgradT = diff(gradT,rho);
% cond = [DgradT(0)== 0, gradT(0)== C];
% gradTSol(rho) = dsolve(eqn,cond)

%% test integral
L =  10; C = 2;
sigma = .05; mu = L/2; K = .14;
% Pdep = K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2);
% x0 = @(rho) C + 0*rho;
x0 = @(rho) K * 1/(sigma*sqrt(pi))*exp(-(1/2)*(rho-mu).^2/sigma.^2) + C;
f1 = @(rho) C*(x0(rho).*cos(2*pi*rho/L));
f2 = @(rho) C*(x0(rho).*cos(2*pi*2*rho/L));
f3 = @(rho) C*(x0(rho).*cos(2*pi*3*rho/L));
f4 = @(rho) C*(x0(rho).*cos(2*pi*4*rho/L));
f5 = @(rho) C*(x0(rho).*cos(2*pi*5*rho/L));
f6 = @(rho) C*(x0(rho).*cos(2*pi*6*rho/L));
f7 = @(rho) C*(x0(rho).*cos(2*pi*7*rho/L));
f8 = @(rho) C*(x0(rho).*cos(2*pi*8*rho/L));
f9 = @(rho) C*(x0(rho).*cos(2*pi*9*rho/L));
f10 = @(rho) C*(x0(rho).*cos(2*pi*10*rho/L));
b0 = 1/L * integral(x0,0,L);
b1 = 2/L * integral(f1,0,L);
b2 = 2/L * integral(f2,0,L);
b3 = 2/L * integral(f3,0,L);
b4 = 2/L * integral(f4,0,L);
b5 = 2/L * integral(f5,0,L);
b6 = 2/L * integral(f6,0,L);
b7 = 2/L * integral(f7,0,L);
b8 = 2/L * integral(f8,0,L);
b9 = 2/L * integral(f9,0,L);
b10 = 2/L * integral(f10,0,L);

T = @(rho,t) exp(-4*0*pi^2*C*t/L^2).*cos(2*pi*0*rho/L)*b0+...
    +exp(-4*1*pi^2*C*t/L^2).*cos(2*pi*1*rho/L)*b1+...
    +exp(-4*2*pi^2*C*t/L^2).*cos(2*pi*2*rho/L)*b2+...
    +exp(-4*3*pi^2*C*t/L^2).*cos(2*pi*3*rho/L)*b3+...
    +exp(-4*4*pi^2*C*t/L^2).*cos(2*pi*4*rho/L)*b4+...
    +exp(-4*5*pi^2*C*t/L^2).*cos(2*pi*5*rho/L)*b5+...
    +exp(-4*6*pi^2*C*t/L^2).*cos(2*pi*6*rho/L)*b6+...
    +exp(-4*7*pi^2*C*t/L^2).*cos(2*pi*7*rho/L)*b7+...
    +exp(-4*8*pi^2*C*t/L^2).*cos(2*pi*8*rho/L)*b8+...
    +exp(-4*9*pi^2*C*t/L^2).*cos(2*pi*9*rho/L)*b9+...
    +exp(-4*10*pi^2*C*t/L^2).*cos(2*pi*10*rho/L)*b10;

fsurf(T,[0 L 0 2])