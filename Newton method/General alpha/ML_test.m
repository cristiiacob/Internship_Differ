function [dTdt] = ML_test(t,T,drho,n,Pdep1,Pdep2,u1,u2,gamma,alpha)

a = 2/3;
s = gamma * a * (alpha + 1) / 2^alpha / (drho)^(alpha+2);

% build the system of ODEs

% Neumann boundary condition
dTdt(1) = s*((T(2) - T(1))^(alpha + 1)) + a*(Pdep1(1)*u1(t) + Pdep2(1)*u2(t));

% Dirichlet boundary condition
dTdt(n+1) = s*0 + a*(Pdep1(n+1)*u1(t) + Pdep2(n+1)*u2(t));

for i = 2:n
    dTdt(i) = s*((T(i+1) - T(i-1))^alpha)*(T(i+1) - 2*T(i) + T(i-1)) + a*(Pdep1(i)*u1(t) + Pdep2(i)*u2(t));
end
dTdt = dTdt';

end


