function [T] = Newton_ML(Told,a,c,drho,dt,Pdep,u)
%   Detailed explanation goes here
s = 1e6 / 4 / (drho)^4; %% is minus ok
q = 2/3 * dt;
n = length(Told) - 1;
F = zeros(n+1,1);
T = Told;
k = 0;
err = 1;
while err > 1e-7 && k <= 100
    F(1) = s*(T(2)^3 - 3*T(1)*T(2)^2 + 3*T(1)^2*T(2) - T(1)^3);
    F(n+1) = 0;
    
    main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
    main(1) = s*(-3*T(2)^2 + 6*T(1)*T(2) - 3*T(1)^2);
    main(n+1) = 0;
    upper(1) = s*(3*T(2)^2 - 6*T(1)*T(2) + 3*T(1)^2);
    
    for i = 2:n 
        F(i) = s*(T(i+1)^3 - 2*T(i)*T(i+1)^2 - T(i-1)*T(i+1)^2 + 4*T(i-1)*T(i)*T(i+1) - (T(i-1)^2)*T(i+1) - 2*(T(i-1)^2)*T(i) + T(i-1)^3);
        main(i) = s*(-2*T(i+1)^2 + 4*T(i-1)*T(i+1) - 2*T(i-1)^2);
        upper(i) = s*(3*T(i+1)^2 - 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) + 4*T(i-1)*T(i) - T(i-1)^2);
        under(i-1) = s*(-T(i+1)^2 + 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) - 4*T(i-1)*T(i) + 3*T(i-1)^2);
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

