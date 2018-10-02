function [dTdt] = nonlin_system(t,T,u,Pdep,b,dx,n)
a = 3/2;
s = 1/a/(2*dx^3);

Aupper = zeros(1,n);
Aunder = zeros(1,n);
Amain = zeros(1,n+1);

for i = 1:n
    Aupper(i) = T(i+1) - 2*T(i);
    Aunder(i) = 2*T(i+1) - T(i);
end
Amain(1) = T(1); 
Aunder(end) = 0;

% A = spdiags([Aunder,Amain,Aupper],[-1 0 1],n+1,n+1);
A = full(gallery('tridiag',n+1,Aunder,Amain,Aupper));
A = s*A;

Bconstr = zeros(n+1,1); Bconstr(1) = -s*b*dx^2;

dTdt = A*T + Bconstr + 1/a*Pdep*u(t);

end