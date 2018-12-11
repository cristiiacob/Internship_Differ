function [dTdt] = ML_general(t,T,u,Pdep,drho,n)
gamma = 1;
alpha = 0;
a = 2/3;
s = gamma * a * (alpha + 1) / 2^alpha / (drho)^(alpha+2);

% build the system of ODEs

% Neumann boundary condition
dTdt(1) = s*((T(2) - T(1))^(alpha + 1)) + a*Pdep(1)*u(t);

% Dirichlet boundary condition
dTdt(n+1) = s*0 + a*Pdep(n+1)*u(t);

for i = 2:n
    dTdt(i) = s*((T(i+1) - T(i-1))^alpha)*(T(i+1) - 2*T(i) + T(i-1)) + a*Pdep(i)*u(t);
end
dTdt = dTdt';

end


