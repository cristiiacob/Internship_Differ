function [F,J] = linear_fsolve(T,Told,C,n,drho,dt,Pdep,u)
gamma = 1;
s = gamma * 1/3 * dt / (drho^3);
% q = 2/3 * dt;
F = zeros(n+1,1);
% F(1) = s*T(2)  - (s)*T(1);
%     F(n+1) =  T(n+1) - C;
%     
%     main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
%     main(1) = -(s);
%     main(n+1) = 1;
%     upper(1) = s;
%     
%     for i = 2:n
%         F(i) = s*T(i+1) + s*T(i-1) - (2*s)*T(i);
%         main(i) = -(2*s);
%         upper(i) = s;
%         under(i-1) = s;
%     end
%     Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
%     Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
%     Under = sparse(2:n+1,1:n,under,n+1,n+1);
%     J = Main + Upper + Under;

F(1) = s*(T(2)^2 + T(1)^2 - 2*T(1)*T(2));
for i = 2:n
    F(i) = s*(T(i+1)^2 - T(i-1)^2 - 2*T(i)*T(i+1) + 2*T(i-1)*T(i));
end
F(n+1) = T(n+1) - C;

% Evaluate the Jacobian if nargout > 1
if nargout > 1
    main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
    %    main(1) = -1;
    main(1) = 2*s*(T(1) - T(2));
    main(n+1) = 1;
    upper(1) = 2*s*(T(2) - T(1));
    
    for i = 2:n
        main(i) = 2*s*(T(i-1) - T(i+1));
        upper(i) = 2*s*(T(i+1) - T(i));
        under(i-1) = 2*s*(T(i) - T(i-1));
    end
    Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
    Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
    Under = sparse(2:n+1,1:n,under,n+1,n+1);
    J = Main + Upper + Under;
    
%     s = 1000000 / 4 / (drho)^4;
%     q = 2/3;
%     n = length(T) - 1;
%     F = zeros(n+1,1);
%     % T = Told;
%     F(1) = s*(T(2)^3 - 3*T(1)*T(2)^2 + 3*T(1)^2*T(2) - T(1)^3);
%     F(n+1) = 0;
%     
%     main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
%     main(1) = s*(-3*T(2)^2 + 6*T(1)*T(2) - 3*T(1)^2);
%     main(n+1) = 0;
%     upper(1) = s*(3*T(2)^2 - 6*T(1)*T(2) + 3*T(1)^2);
%     
%     for i = 2:n
%         F(i) = s*(T(i+1)^3 - 2*T(i)*T(i+1)^2 - T(i-1)*T(i+1)^2 + 4*T(i-1)*T(i)*T(i+1) - (T(i-1)^2)*T(i+1) - 2*(T(i-1)^2)*T(i) + T(i-1)^3);
%         main(i) = s*(-2*T(i+1)^2 + 4*T(i-1)*T(i+1) - 2*T(i-1)^2);
%         upper(i) = s*(3*T(i+1)^2 - 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) + 4*T(i-1)*T(i) - T(i-1)^2);
%         under(i-1) = s*(-T(i+1)^2 + 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) - 4*T(i-1)*T(i) + 3*T(i-1)^2);
%     end
%     Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
%     Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
%     Under = sparse(2:n+1,1:n,under,n+1,n+1);
%     J = Main + Upper + Under;

%     s = 1000000 * dt / 4 / (drho)^4;
%     q = 2/3;
%     n = length(T) - 1;
%     F = zeros(n+1,1);
%     % T = Told;
%     F(1) = s*(T(2)^3 - 3*T(1)*T(2)^2 + 3*T(1)^2*T(2) - T(1)^3) - T(1) + Told(1);
%     F(n+1) = T(n+1) - C;
%     
%     main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
%     main(1) = s*(-3*T(2)^2 + 6*T(1)*T(2) - 3*T(1)^2) - 1;
%     main(n+1) = 1;
%     upper(1) = s*(3*T(2)^2 - 6*T(1)*T(2) + 3*T(1)^2);
%     
%     for i = 2:n
%         F(i) = s*(T(i+1)^3 - 2*T(i)*T(i+1)^2 - T(i-1)*T(i+1)^2 + 4*T(i-1)*T(i)*T(i+1) - (T(i-1)^2)*T(i+1) - 2*(T(i-1)^2)*T(i) + T(i-1)^3) - T(i) + Told(i);
%         main(i) = s*(-2*T(i+1)^2 + 4*T(i-1)*T(i+1) - 2*T(i-1)^2) - 1;
%         upper(i) = s*(3*T(i+1)^2 - 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) + 4*T(i-1)*T(i) - T(i-1)^2);
%         under(i-1) = s*(-T(i+1)^2 + 4*T(i)*T(i+1) - 2*T(i-1)*T(i+1) - 4*T(i-1)*T(i) + 3*T(i-1)^2);
%     end
%     Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
%     Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
%     Under = sparse(2:n+1,1:n,under,n+1,n+1);
%     J = Main + Upper + Under;
end

