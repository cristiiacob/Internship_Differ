function [dTdt] = nonlin_sys_3rd(t,T,u,Pdep,b,dx,n)

s = 1000000 / 4 / (dx)^4; %% is minus ok
q = 2/3;
    dTdt(1) = s*(T(2)^3 - 3*T(1)*T(2)^2 + 3*T(1)^2*T(2) - T(1)^3) + q*Pdep(1)*u(t);
    dTdt(n+1) = s*0 + q*Pdep(n+1)*u(t);
    for i = 2:n 
       dTdt(i) = s*(T(i+1)^3 - 2*T(i)*T(i+1)^2 - T(i-1)*T(i+1)^2 + 4*T(i-1)*T(i)*T(i+1) - ...
             (T(i-1)^2)*T(i+1) - 2*(T(i-1)^2)*T(i) + T(i-1)^3) + q*Pdep(i)*u(t);
    end
    dTdt = dTdt';

end