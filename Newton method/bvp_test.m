function bvp_test
tic
    rho_start = 0; rho_end = 1000; n = 30000;
    solinit = bvpinit(linspace(rho_start,rho_end,n+1),[2 .15]);
    options = bvpset('NMax',30000);
    sol = bvp5c(@bvp_ode,@bvp_bc,solinit,options);
    rho_int = linspace(rho_start,rho_end,100);
    Srho_int = deval(sol,rho_int);
    plot(rho_int,Srho_int(1,:));
    toc
end

%------------------------------------------
function dTdrho = bvp_ode(rho,T)
P = .000001; alpha = 1/2;
dTdrho = [T(2) -(P/(alpha+1))*1/((T(2)^alpha))];
end

function res = bvp_bc(Ta,Tb)
C = 2;
res = [Ta(2) Tb(1)-C];
end