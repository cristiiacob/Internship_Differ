function [dTdt] = ML_1(t,T,u,Pdep,drho,n)

s = 1 / 4 / (drho)^2;
q = 2 / (drho)^3;

% build the system of ODEs

% Neumann boundary condition
dTdt(1) = s*(1/T(1)^2)*(T(2)^2-2*T(1)*T(2)+T(1)^2) - q*(1/T(1))*(T(2)^2-2*T(1)*T(2)+T(1)^2) + Pdep(1)*u(t);

% Dirichlet boundary condition
dTdt(n+1) = s*0 + Pdep(n+1)*u(t);

for i = 2:n
    dTdt(i) = s*(1/T(i)^2)*(T(i+1)^2-2*T(i-1)*T(i+1)+T(i-1)^2) ...
        -q*(1/T(i))*(T(i+1)^2-2*T(i)*T(i+1)+2*T(i-1)*T(i)-T(i-1)^2) + Pdep(i)*u(t);
end
dTdt = dTdt';

end


