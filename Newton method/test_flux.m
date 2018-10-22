gradT = -(0:.0001:.02);
chi = 1/2;
alpha = 1/2;

q = (-chi * gradT) .^ (alpha+1);

plot(-gradT,q);
%%
syms gradT_sym
q_sym = (-chi * gradT_sym) ^ (alpha+1);

pol_taylor = taylor(q_sym,gradT_sym,'ExpansionPoint',-0.01,'Order',4)
pol_coeffs = collect(pol_taylor,gradT_sym)
q_fitted = subs(pol_coeffs,gradT_sym,gradT);
plot(-gradT,[q; q_fitted]);



%% test initial temperature profile - not finding solutions
% syms gradT(rho) P alpha C
% eqn = diff(gradT,rho,2) == -P / (alpha+1) * 1/(gradT ^ alpha);
% DgradT = diff(gradT,rho);
% cond = [DgradT(0)== 0, gradT(0)== C];
% gradTSol(rho) = dsolve(eqn,cond)