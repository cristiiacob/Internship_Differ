function [T] = linear_fct(Told,a,c,drho,dt,Pdep,u)
%   Detailed explanation goes here
s = 2/3 * dt / (drho)^2;
q = 2/3 * dt;
n = length(Told) - 1;
F = zeros(n+1,1);
T = Told;
k = 0;
err = 1;
while err > 1e-7 && k <= 100
    F(1) = Told(1) + s*T(2)  - (s+1)*T(1)  + q*Pdep(1)*u;
    F(n+1) =  T(n+1) - c + q*Pdep(n+1)*u;
    
    main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
    main(1) = -(s + 1);
    main(n+1) = 1;
    upper(1) = s;
    
    for i = 2:n
        F(i) = Told(i) + s*T(i+1) + s*T(i-1) - (2*s + 1)*T(i) + q*Pdep(i)*u;
        main(i) = -(2*s + 1);
        upper(i) = s;
        under(i-1) = s;
    end
    Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
    Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
    Under = sparse(2:n+1,1:n,under,n+1,n+1);
    J = Main + Upper + Under;
    
    % Newton step
    Tprevious = T;
    step = -J\F;
    T = T + step;
    k = k+1;
    err = norm(T - Tprevious);
end
end

