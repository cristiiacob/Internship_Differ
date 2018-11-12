function [dTdt] = lin_eq(t,T,u,Pdep,b,dx,n)
    gamma = 1/2;
    a = 3/2; 
    s = gamma * 1/a/(dx^2);
    dTdt (1) = s*(T(2) - T(1) - 2*b*dx) + 1/a*Pdep(1)*u(t); % rmove 2 in fron of t(2) and T(1)?
    dTdt(n+1) = s*0 + 1/a*Pdep(n+1)*u(t);
    for i = 2:n
       dTdt(i) = s*(T(i+1) - 2*T(i) + T(i-1)) + 1/a*Pdep(i)*u(t); 
    end
    dTdt = dTdt';
end