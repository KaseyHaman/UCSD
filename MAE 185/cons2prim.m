function [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv)
%Converts conservative variables to primitive variables

%Solve for variables
rho = squeeze(U(1, :, :));
u = squeeze(U(2, :, :))./rho;
v = squeeze(U(3, :, :))./rho;
Et = squeeze(U(4, :, :));
e = (Et ./ rho) - ((u.^2 +v.^2).*0.5);
T = e ./ cv;
p = rho.*R.*T;


end