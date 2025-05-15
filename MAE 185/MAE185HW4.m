%% 4.1
clc;clear all;
%Load data
load('supersonicJetLES_xyPlane')

%Define variables
R = 287; %J/kg K
cp = 1005; %J/kg K
cv = 718; %J / kg K

%Prim2cons 
U = prim2cons(rho, u, v, T, cv);

%Cons2prim
[rho2, u2, v2, T2, p2, e2, Et2] = cons2prim(U, R, cv);

%Sutherland
[mu] = sutherland(T);

%Check step account for floating point errors
error = 1e-5;
if all(abs(u2(:) - u(:)) < error) && all(abs(T2(:) - T(:)) < error) && all(abs(v2(:) - v(:)) < error)
    disp('Variables match.')
else 
    disp('Variables do not match.')
end

%Define pu and pv seperate
pu = squeeze(U(2, :, :));
pv = squeeze(U(3, :, :));

%List out all variables for subplots
variables = {rho, pu, pv, T, p2, e2, Et2, mu};
titles = {'\rho', '\rhou', '\rhov', 'T', 'p', 'e', 'E_t', '\mu'};
units = {'[kg/m^3]', '[kg/(m^2·s)]', '[kg(m^2·s)]', '[K]', '[Pa]', '[J/kg]', '[J/m^3]', '[N·s/m^2]'};

%Create a for loop to solve the system
for i = 1:length(variables)
    subplot(4, 2, i)
    pcolor(xx, yy, variables{i}) 
    shading interp
    axis equal tight
    title(titles{i})
    xlabel('x')
    ylabel('y')
    cb = colorbar;
    ylabel(cb, [titles{i} ' ' units{i}])
end


%% 4.2
clc;clear all;
%Load data
load('supersonicJetLES_xyPlane')

%Define Step Sizes
dx = xx(2, 1)-xx(1, 1);
dy = yy(1, 2)-yy(1, 1);

%Define variables
R = 287; %J/kg K
cp = 1005; %J/kg K
cv = 718; %J / kg K

%Prim2cons 
U = prim2cons(rho, u, v, T, cv);

%Cons2prim
[rho, u, v, T, p, e, Et] = cons2prim(U, R, cv);

%Sutherland
[mu] = sutherland(T);

%Forward differences in x and central differences in y
Txx1 = 2.*mu.*(ddx_fwd(u, dx)-1/3.*(ddx_fwd(u, dx) + ddy_central(v, dy)));
Tyy1 = 2.*mu.*(ddy_central(v, dy)-1/3.*(ddx_fwd(u, dx) + ddy_central(v, dy)));
Txy1 = mu.*(ddy_central(u, dy)+ddx_fwd(v, dx));

%Central differences in x and backward differences in y
Txx2 = 2.*mu.*(ddx_central(u, dx)-1/3.*(ddx_central(u, dx) + ddy_bwd(v, dy)));
Tyy2 = 2.*mu.*(ddy_bwd(v, dy)-1/3.*(ddx_central(u, dx) + ddy_bwd(v, dy)));
Txy2 = mu.*(ddy_bwd(u, dy)+ddx_central(v, dx));

%Sort into variable
shearstress = {Txx1, Tyy1, Txy1, Txx2, Tyy2, Txy2};
titles = {'\tau_{xx} #1', '\tau_{yy} #1', '\tau_{xy} #1', '\tau_{xx} #2', '\tau_{yy} #2', '\tau_{xy} #2'};
units = '[N / m^2]';

for i = 1:length(shearstress)
    subplot(2, 3, i)
    pcolor(xx, yy, shearstress{i}) 
    shading interp
    axis equal tight
    clim([-0.5 0.5])  
    title(titles{i})
    xlabel('x')
    ylabel('y')
    colorbar;
    ylabel(colorbar, [titles{i} ' ' units])
end



%% 4.3
clc;clear all;
%Load data
load('gotritons_uv')

%Define Step Sizes
dx = xx(2, 1)-xx(1, 1);
dy = yy(1, 2)-yy(1, 1);

%Solve for first time step
cx = max(max((u)));
cy = max(max((v)));
SF = 1.1;
dt_limit = (dx*dy) / (cx*dy + cy*dx);
dt = dt_limit / SF;

%Define variables and paramaters
d = 0.0005;
bl = -0.5;
bu = 0.5;
t = 0;
t_final = 2;

%Maintain boundary conditions
    u(:, 1) = 0;
    u(:, end) = 0;
    v(:, 1) = 0;
    v(:, end) = 0;
%Run Models
while t < t_final

    pcolor(xx, yy, u);
    shading interp;
    axis equal tight;
    title(['u @ t=' num2str(t)])
    xlabel('x');
    ylabel('y');
    colorbar;
    clim([0 2]);
    drawnow;

    %Solve for predictor time step using MacCormack method
    u_pred = u + dt * (-u.*ddx_fwd_periodic(u, dx) - v.*ddy_fwd(u, dy) + d.*(d2dx2periodic(u, dx) + d2dy2(u, dy))); %Set new u_predictor value

    v_pred = v + dt * (-u.*ddx_fwd_periodic(v, dx) - v.*ddy_fwd(v, dy) + d.*(d2dx2periodic(v, dx) + d2dy2(v, dy))); %Set new v_predictor value

    %Maintain boundary conditions
    u_pred(:, 1) = 0;
    u_pred(:, end) = 0;
    v_pred(:, 1) = 0;
    v_pred(:, end) = 0;
    
    %Solve for next time step using MacCormack method
    u = 1/2.*(u + u_pred) + 1/2 * dt * (-u_pred.*ddx_bwd_periodic(u_pred, dx) - v_pred.*ddy_bwd(u_pred, dy) + d.*(d2dx2periodic(u_pred, dx) + d2dy2(u_pred, dy)));

    v = 1/2.*(v + v_pred) + 1/2 * dt * (-u_pred.*ddx_bwd_periodic(v_pred, dx) - v_pred.*ddy_bwd(v_pred, dy) + d.*(d2dx2periodic(v_pred, dx) + d2dy2(v_pred, dy)));

    t = t + dt; %increase time

    %Maintain boundary conditions
    u(:, 1) = 0;
    u(:, end) = 0;
    v(:, 1) = 0;
    v(:, end) = 0;

    %Set new stability time step before next iteration
    cx = max(max((u)));
    cy = max(max((v)));
    SF = 1.1;
    dt_limit = (dx*dy) / (cx*dy + cy*dx);
    dt = dt_limit / SF;
   
end





