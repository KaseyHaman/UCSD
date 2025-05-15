%Write a MATLAB code that solves the 2-D compressible Navier-Stokes equations with closure for calorically perfect gas using MacCormack's method. 
%To demonstrate your code, solve the problem of the flow over a flat plate with a sharp leading-edge at Mach 4. 

%% Part 1

clc
clear
adiabatic = 0;

%Physical Parameters
M = 4;
L = 1e-5; %m
H = 8e-6; %m
cv = 718; %J/kg/K
cp = 1005; %J/kg/K
R = cp-cv;
p_atm = 101300; %Pa
T_atm = 288.15; %K
gamma = 1.4;
a = sqrt(gamma*R*T_atm); %Speed of sound of air at atm conditions (m/s)



%Specify number of grid points
nx = 75;
ny = 80;

%Initialize Gridpoints
[xx,yy] = ndgrid(linspace(0,L,nx),linspace(0,H,ny));

%Find dx,dy
dx = xx(2,1)-xx(1,1);
dy = yy(1,2)-yy(1,1);

%Use given timestep
dt = 2.35e-11;

%Initialize primitive variables

%Set primitive variable fields with IC's
p = ones(size(xx))*p_atm;
T = ones(size(xx))*T_atm;
rho = p./(R*T);
u = ones(size(xx))*M*a;
v = zeros(size(xx));


%Initialize conservative vars U vector
U = prim2cons(rho,u,v,T,cv);

%Initialize conservative variables
U_pred = zeros(size(U));
E = zeros(size(U));
F = zeros(size(U));

%Initialize Convergence Matrix of u
u_conv = zeros(1,1500);

%Set BC's
U = enforcebcs(U, adiabatic);

%Get Updated Primitive Variables from U
[rho,u,v,T,p,e,~] = cons2prim(U,R,cv);

%Plot IC for rho,u,v,e,p,T
figure(1)
tiledlayout(2,3)
nexttile
pcolor(xx,yy,rho),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'\rho [kg/m^3]')
title('\rho')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,u),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'u [m/s]')
title('u')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,v),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'v [m/s]')
title('v')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,e),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'e [J]')
title('e')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,p),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'p [Pa]')
title('p')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,T),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'T [K]')
title('T')
xlabel('x')
ylabel('y')

%% 

%Start time loop for 1500 time steps
for n = 1:1500

u_n = u;

%Predictor

%Compute flux vectors for predictor step
E = xfluxvec_if_fwd(U,dx,dy);
F = yfluxvec_if_fwd(U,dx,dy);

%Compute U_pred using fwd differences
U_pred = U - dt*fd_for_3d_matrix(E,dx,dy,'ddx_fwd') - dt*fd_for_3d_matrix(F,dx,dy,'ddy_fwd');

%Enforce BC's on predictor U
U_pred = enforcebcs(U_pred, adiabatic);

%Compute flux vectors for corrector step
E = xfluxvec_if_bwd(U_pred,dx,dy);
F = yfluxvec_if_bwd(U_pred,dx,dy);

%Corrector Step
U = (1/2)*(U_pred + U) + (1/2)*(- dt*fd_for_3d_matrix(E,dx,dy,'ddx_bwd') - dt*fd_for_3d_matrix(F,dx,dy,'ddy_bwd'));

%Enforce BC's on corrector U
U = enforcebcs(U, adiabatic);

%Get Updated Primitive Variables from U
[rho,u,v,T,p,e,~] = cons2prim(U,R,cv);


%Plots
if mod(n,50) == 0
%Plot IC
figure(1)
tiledlayout(2,3)
nexttile
pcolor(xx,yy,rho),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'\rho [kg/m^3]')
title('\rho')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,u),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'u [m/s]')
title('u')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,v),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'v [m/s]')
title('v')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,e),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'e [J]')
title('e')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,p),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'p [Pa]')
title('p')
xlabel('x')
ylabel('y')

nexttile
pcolor(xx,yy,T),shading interp, axis equal tight
cb = colorbar; ylabel(cb,'T [K]')
title('T')
xlabel('x')
ylabel('y')

drawnow;


end

disp(max(max(u_n-u)));
u_conv(n) = max(max(u_n-u));
figure(2)
plot(1:1500,u_conv)
end


%Plot visuals for part 2.3
x_locations = round([0.25; 0.5; 0.75] .* nx);

figure(3)
for i = 1:length(x_locations)
subplot(2, 1, 1)
    plot(T(x_locations(i), :) / T_atm, yy(x_locations(i), :));
hold on;
title('Adiabatic Normalized Wall Temperature at Key Points')
xlabel('Normalized Temperature')
ylabel('y')

subplot(2, 1, 2)
 plot(p(x_locations(i), :) / p_atm, yy(x_locations(i), :));
hold on;
title('Normalized Wall Pressure at Key Points')
xlabel('Normalized Pressure')
ylabel('y')

end
legend('x/L = 0.25', 'x/L = 0.5', 'x/L=0.75')

%Plot visuals for part 2.3.2
figure(4)
plot(xx(:, 1), T(:, 1), 'k', 'linewidth', 1.5);
title('Wall Temperature as a Function of x')
xlabel('x')
ylabel('Wall Temperature')


%% Part 2.1

%Define constants
beta = 0.8;
kappa = 10;

%Find the magnitude of the gradient of density
drhodx = ddx_central(rho, dx);
drhody = ddy_central(rho, dy);
mag_grad_rho = sqrt(drhodx.^2 + drhody.^2);

%Compute better results from scaled quantity
S = beta * exp(-kappa * mag_grad_rho / max(mag_grad_rho(:)));

%Plot our results
figure
pcolor(xx, yy, S), shading interp, axis equal tight
colormap(gray)
colorbar
title('Numerical Schlieren')
xlabel('x'); ylabel('y')
hold on;

%% Part 2.2
theta = asin(1/M);
x_shock = linspace(0, L, 100);  
y_shock = tan(theta) * x_shock;  
plot(x_shock, y_shock, 'r', 'LineWidth', 1.5);


%% Part 2.3 (Redo part 1 with adiabatic)
clc
clear

adiabatic = 1;

%Physical Parameters
M = 4;
L = 1e-5; %m
H = 8e-6; %m
cv = 718; %J/kg/K
cp = 1005; %J/kg/K
R = cp-cv;
p_atm = 101300; %Pa
T_atm = 288.15; %K
gamma = 1.4;
a = sqrt(gamma*R*T_atm); %Speed of sound of air at atm conditions (m/s)



%Specify number of grid points
nx = 75;
ny = 80;

%Initialize Gridpoints
[xx,yy] = ndgrid(linspace(0,L,nx),linspace(0,H,ny));

%Find dx,dy
dx = xx(2,1)-xx(1,1);
dy = yy(1,2)-yy(1,1);

%Use given timestep
dt = 2.35e-11;

%Initialize primitive variables

%Set primitive variable fields with IC's
p = ones(size(xx))*p_atm;
T = ones(size(xx))*T_atm;
rho = p./(R*T);
u = ones(size(xx))*M*a;
v = zeros(size(xx));


%Initialize conservative vars U vector
U = prim2cons(rho,u,v,T,cv);

%Initialize conservative variables
U_pred = zeros(size(U));
E = zeros(size(U));
F = zeros(size(U));

%Initialize Convergence Matrix of u
u_conv = zeros(1,1500);

%Set BC's
U = enforcebcs(U, adiabatic);

%Get Updated Primitive Variables from U
[rho,u,v,T,p,e,~] = cons2prim(U,R,cv);


%Start time loop for 1500 time steps
for n = 1:1500

u_n = u;

%Predictor

%Compute flux vectors for predictor step
E = xfluxvec_if_fwd(U,dx,dy);
F = yfluxvec_if_fwd(U,dx,dy);

%Compute U_pred using fwd differences
U_pred = U - dt*fd_for_3d_matrix(E,dx,dy,'ddx_fwd') - dt*fd_for_3d_matrix(F,dx,dy,'ddy_fwd');

%Enforce BC's on predictor U
U_pred = enforcebcs(U_pred, adiabatic);

%Compute flux vectors for corrector step
E = xfluxvec_if_bwd(U_pred,dx,dy);
F = yfluxvec_if_bwd(U_pred,dx,dy);

%Corrector Step
U = (1/2)*(U_pred + U) + (1/2)*(- dt*fd_for_3d_matrix(E,dx,dy,'ddx_bwd') - dt*fd_for_3d_matrix(F,dx,dy,'ddy_bwd'));

%Enforce BC's on corrector U
U = enforcebcs(U, adiabatic);

%Get Updated Primitive Variables from U
[rho,u,v,T,p,e,~] = cons2prim(U,R,cv);


end


%% Plot specific visuals for part 2.3.1
x_locations = round([0.25; 0.5; 0.75] .* nx);

figure(1)
for i = 1:length(x_locations)
subplot(2, 1, 1)
    plot(T(x_locations(i), :) / T_atm, yy(x_locations(i), :));
hold on;
title('Adiabatic Normalized Wall Temperature at Key Points')
xlabel('Normalized Temperature')
ylabel('y')

subplot(2, 1, 2)
 plot(p(x_locations(i), :) / p_atm, yy(x_locations(i), :));
hold on;
title('Normalized Wall Pressure at Key Points')
xlabel('Normalized Pressure')
ylabel('y')

end
legend('x/L = 0.25', 'x/L = 0.5', 'x/L=0.75')

%% Plot visuals for part 2.3.2
figure(2)
plot(xx(:, 1), T(:, 1), 'k', 'linewidth', 1.5);
title('Wall Temperature as a Function of x')
xlabel('x')
ylabel('Wall Temperature')


%% Functions
function U_bc = enforcebcs(U, adiabatic)

cv = 718; %J/kg/K
cp = 1005; %J/kg/K
R = cp-cv;
p_atm = 101300; %Pa
T_atm = 288.15; %K
gamma = 1.4;
a = sqrt(gamma*R*T_atm); %Speed of sound of air at atm conditions (m/s)
M = 4;
u_inf = M*a;

[~,u,v,T,p,~,~] = cons2prim(U,R,cv);

%Enforce BC on u
u(1,:) = u_inf;
u(:,1) = 0;
u(end,:) = 2*u(end-1,:)-u(end-2,:);
u(:,end) = u_inf;

%Enforce BC on v
v(1,:) = 0;
v(:,1) = 0;
v(end,:) = 2*v(end-1,:)-v(end-2,:);
v(:,end) = 0;

%Enforce BC on p
p(1,:) = p_atm;
p(:,end) = p_atm;
p(end,:) = 2*p(end-1,:)-p(end-2,:);
p(:,1) = 2*p(:,2)-p(:,3);

%Enforce BC on T
T(1,:) = T_atm;
T(:,end) = T_atm;
T(end,:) = 2*T(end-1,:)-T(end-2,:);
if adiabatic == 0
    T(:,1) = T_atm;
else
    T(:, 1) = T(:, 2);
end

%Find rho field
rho = p./(R*T);

%Get new U vector
U_bc = prim2cons(rho,u,v,T,cv);

end