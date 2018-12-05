function [T] = Newton_test(Told,c,drho,dt,Pdep1,Pdep2,u1,u2,gamma,alpha)
a = 2/3;
s = gamma * dt * a * (alpha + 1) / 2^alpha / (drho)^(alpha+2);
q = a * dt;
n = length(Told) - 1;
F = zeros(n+1,1);
T = Told;
k = 0;
err = 1;

% main loop for the Newton iteration algorithm
while err > 1e-7 && k <= 100
    % build the system of equations F and the Jacobain J
    
    % Neumann boundary condition
    F(1) = s*((T(2) - T(1))^(alpha + 1)) - T(1) + Told(1) + q*(Pdep1(1)*u1 + Pdep2(1)*u2);
    % Dirichlet boundary condition
    F(n+1) =  T(n+1) - c + q*(Pdep1(n+1)*u1 + Pdep2(n+1)*u2);
    
    % main diagonal, upper diagonal adn under diagonal
    main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
    main(1) = -s*(alpha + 1)*((T(2) - T(1))^alpha) - 1;
    main(n+1) = 1;
    upper(1) = s*(alpha + 1)*((T(2) - T(1))^alpha);
    
    if alpha == 0
        for i = 2:n
            F(i) = s*T(i+1) + s*T(i-1) - (2*s + 1)*T(i) + Told(i) + q*(Pdep1(i)*u1 + Pdep2(i)*u2);
            main(i) = -2*s - 1;
            upper(i) = s;
            under(i-1) = s;
        end
    else
        for i = 2:n
            F(i) = s*((T(i+1) - T(i-1))^alpha)*(T(i+1) - 2*T(i) + T(i-1)) - T(i) + Told(i) + q*(Pdep1(i)*u1 + Pdep2(i)*u2);
            main(i) = -2*s*((T(i+1) - T(i-1))^alpha) - 1;
            upper(i) = s*((T(i+1) - T(i-1))^(alpha - 1))*((alpha + 1)*T(i+1) - 2*alpha*T(i) + (alpha - 1)*T(i-1));
            under(i-1) = s*((T(i+1) - T(i-1))^(alpha - 1))*((-alpha + 1)*T(i+1) + 2*alpha*T(i) - (alpha + 1)*T(i-1));
        end
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

