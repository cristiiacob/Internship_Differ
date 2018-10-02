function [dTdt] = nonlin_eq(t,T,u,Pdep,b,dx,n)
a = 3/2;
s = 1/a/(2*dx^3);
%     s = 1;
%     dTdt (1) = s*(4*b*dx*T(2) - 4*b*dx*T(1) - 4*b^2*dx^2) + 1/a*Pdep(1)*u(t);
dTdt (1) = s*(T(2)^2 + T(1)^2 - T(1)*T(2) - b^2*dx^2) + 1/a*Pdep(1)*u(t);
dTdt(n+1) = s*0 + 1/a*Pdep(n+1)*u(t);
% dTdt(n+1) = 1;
%     dTdt(n+1) = s*4*(T(n)*b*dx + b^2*dx^2 + T(n+1)*b*dx + 4*T(n)*T(n+1));
%     if a dT/dr for ro_end is also a boundary condition
for i = 2:n
    dTdt(i) = s*(T(i+1)^2 - T(i-1)^2 - 2*T(i)*T(i+1) + 2*T(i-1)*T(i)) + 1/a*Pdep(i)*u(t);
end
dTdt = dTdt';
end