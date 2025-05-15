function [E,mu] = xfluxvec_if_fwd(U,dx,dy)
E = zeros(size(U));


cv = 718; %J/kg/K
cp = 1005; %J/kg/K
R = cp-cv;

[~,u,v,T,p,~,Et] = cons2prim(U,R,cv);

%Update mu and k fields
mu = sutherland(T);
k = mu*cp/0.71;

%Compute Shear Stresses
tau_xx = 2*mu.*(ddx_bwd(u,dx)-(1/3)*(ddx_bwd(u,dx)+ddy_central(v,dy)));
tau_xy = mu.*(ddy_central(u,dy)+ddx_bwd(v,dx));

%Compute Heat flux
qx = -k.*ddx_bwd(T,dx);

%Compute X Flux Vector
E(1,:,:) = squeeze(U(2,:,:));

E(2,:,:) = squeeze(U(2,:,:)).*u + p - tau_xx;

E(3,:,:) = squeeze(U(2,:,:)).*v - tau_xy;

E(4,:,:) = (Et + p).*u - u.*tau_xx - v.*tau_xy + qx;

end

