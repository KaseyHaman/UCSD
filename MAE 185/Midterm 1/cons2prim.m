function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
%Function to convert conservative variables to primative variables
%Takes in array U with primitive variables as well as scalar values for R
%and cv
%Outputs all conservative variable fields of size nx by ny

%Convert U array into usable grid variables
rho = squeeze(U(1,:,:));
rho_u = squeeze(U(2,:,:));
rho_v = squeeze(U(3,:,:));
Et = squeeze(U(4,:,:));

%Solve for primitive variables
u = rho_u./rho;
v = rho_v./rho;
e = Et./rho - (u.^2+v.^2)/2;
T = e/cv;
p = rho.*T*R;

end

