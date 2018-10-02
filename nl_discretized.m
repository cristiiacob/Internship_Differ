function [F,J] = nl_discretized(T,Told,C,n,drho,dt)

s = 2/3 * dt / (drho^3);

F = zeros(n+1,1);
% F(1) = Told(1) - T(1);
F(1) = Told(1) + s*T(2)^2 + s*T(1)^2 - 2*s*T(1)*T(2) - T(1);
for i = 2:n
    F(i) = Told(i) + s*(T(i+1))^2-s*(T(i-1))^2 - T(i)*(2*s*T(i+1) - 2*s*T(i-1) + 1);
end
F(n+1) = T(n+1) - C;

% Evaluate the Jacobian if nargout > 1
if nargout > 1
   main = zeros(n+1,1); upper = zeros(n,1); under = zeros(n,1);
%    main(1) = -1; 
   main(1) = 2*s*(T(1) - T(2)) - 1;
   main(n+1) = 1;
   upper(1) = 2*s*(T(2) - T(1));
   
   for i = 2:n
      main(i) = 2*s*(T(i-1) - T(i+1)) - 1;
      upper(i) = 2*s*(T(i+1) - T(i));
      under(i-1) = 2*s*(T(i) - T(i-1));
   end
   Main = sparse(1:n+1,1:n+1,main,n+1,n+1);
   Upper = sparse(1:n,2:n+1,upper,n+1,n+1);
   Under = sparse(2:n+1,1:n,under,n+1,n+1);
   J = Main + Upper + Under;
end

end

