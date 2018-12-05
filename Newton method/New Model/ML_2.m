function [dTdt] = ML_2(t,T,u,Pdep,drho,n)
gamma = 1;
a = 2/3;
s = gamma * a * 3 / 4 / (drho)^4;
q = gamma * a * 1 / (drho)^3;

% build the system of ODEs

% Neumann boundary condition
dTdt(1) = s*(1/T(1)^2)*(T(2)^3-3*T(1)*T(2)^2+3*T(1)^2*T(2)-T(1)^3) - q*(1/T(1)^3)*(T(2)^3-3*T(1)*T(2)^2+3*T(1)^2*T(2)-T(1)^3) + a*Pdep(1)*u(t);

% Dirichlet boundary condition
dTdt(n+1) = s*0 + Pdep(n+1)*u(t);

for i = 2:n
    dTdt(i) = s*(1/T(i)^2)*(T(i+1)^3-2*T(i)*T(i+1)^2-T(i-1)*T(i+1)^2+4*T(i-1)*T(i)*T(i+1)-T(i-1)^2*T(i+1)-2*T(i-1)^2*T(i)+T(i-1)^3) ...
        -q*(1/T(i)^3)*(T(i+1)^3-T(i-1)^3-3*T(i-1)*T(i+1)^2+3*T(i-1)^2*T(i+1)) + a*Pdep(i)*u(t);
end
dTdt = dTdt';

end


