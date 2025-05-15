function F = yfluxvec_if_fwd(U,dx,dy)
F = zeros(size(U));


cv = 718; %J/kg/K
cp = 1005; %J/kg/K
R = cp-cv;

[~,u,v,T,p,~,Et] = cons2prim(U,R,cv);

%Update mu and k fields
mu = sutherland(T);
k = mu*cp/0.71;

%Compute Shear Stresses
tau_yy = 2*mu.*(ddy_bwd(v,dy)-(1/3)*(ddx_central(u,dx)+ddy_bwd(v,dy)));
tau_xy = mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));

%Compute Heat flux
qy = -k.*ddy_bwd(T,dy);

%Compute y Flux Vector
F(1,:,:) = squeeze(U(3,:,:));

F(2,:,:) = squeeze(U(3,:,:)).*u - tau_xy;

F(3,:,:) = squeeze(U(3,:,:)).*v + p - tau_yy;

F(4,:,:) = (Et + p).*v - v.*tau_yy - u.*tau_xy + qy;

end



