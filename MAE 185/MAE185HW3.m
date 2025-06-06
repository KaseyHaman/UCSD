%% 3.1
clc;clear all;
%Load data
load('gotritons','T','xx','yy')

%Define Step Sizes
dx = xx(2, 1)-xx(1, 1);
dy = yy(1, 2)-yy(1, 1);

%Define variables
alpha = 2;
SF = 1.1;
dt_limit = 1/4 / (alpha * (1/dx^2 + 1/dy^2));
dt = dt_limit / SF;
t_final = 0.001;
t = 0;

%Run Model
while t < t_final

  pcolor(xx, yy, T);
    shading interp;
    axis equal tight;
title(['T @ t=' num2str(t)])
    xlabel('x');
    ylabel('y');

    colorbar;
    clim([0 1]);
    drawnow;

    T = T + dt * alpha * (d2dx2periodic(T, dx) + d2dy2periodic(T, dy)); %Set new T

    t = t + dt; %increase time

   
end

%% 3.2
clc;clear all;
%Load data
load('gotritons','T','xx','yy')

%Define Step Sizes
dx = xx(2, 1)-xx(1, 1);
dy = yy(1, 2)-yy(1, 1);

%Define variables
cx = 1;
cy = 1;
SF = 1.1;
dt_limit = 1 / (cx/dx + cy/dy);
dt = dt_limit / SF;

t_final_bwd = 2;
t_final_central = 0.25;
t_bwd = 0;
t_central = 0;

%Run Models
figure(1)
while t_bwd < t_final_bwd

  pcolor(xx, yy, T);
    shading interp;
    axis equal tight;
title(['T @ t=' num2str(t_bwd)])
    xlabel('x');
    ylabel('y');

    colorbar;
    clim([0 1]);
    drawnow;

    T = T + dt * (-ddx_bwd_periodic(T, dx) - ddy_bwd_periodic(T, dy)); %Set new T

    t_bwd = t_bwd + dt; %increase time

   
end

%Clear T and reload data
clear T;
load('gotritons','T','xx','yy')
figure(2);

while t_central < t_final_central

  pcolor(xx, yy, T);
    shading interp;
    axis equal tight;
title(['T @ t=' num2str(t_central)])
    xlabel('x');
    ylabel('y');
    colorbar;
    clim([0 1]);
    drawnow;

    T = T + dt * (-ddx_central_periodic(T, dx) - ddy_central_periodic(T, dy)); %Set new T

    t_central = t_central + dt; %increase time

   
end


%% 3.3
clc;clear all;
%Load data
load('gotritons','T','xx','yy')

%Define Step Sizes
dx = xx(2, 1)-xx(1, 1);
dy = yy(1, 2)-yy(1, 1);

%Define variables
cx = 1;
cy = 1;
SF = 1.1;
dt_limit = 1 / (cx/dx + cy/dy);
dt = dt_limit / SF;
T_pred = 0;
t = 0;
t_final = 2;


%Run Models
while t < t_final

  pcolor(xx, yy, T);
    shading interp;
    axis equal tight;
title(['T @ t=' num2str(t)])
    xlabel('x');
    ylabel('y');

    colorbar;
    clim([0 1]);
    drawnow;

    T_pred = T + dt * (-ddx_fwd_periodic(T, dx) - ddy_fwd_periodic(T, dy)); %Set new T_predictor value

    T = (1/2 * (T + T_pred)) - (1/2 * dt * (ddx_bwd_periodic(T_pred, dx) + ddy_bwd_periodic(T_pred, dy)));

    t = t + dt; %increase time

   
end

%% 3.1 Prof
clc;clear all
load('gotritons.mat','T','xx','yy');
% Solver parameters
alpha = 2; % thermal diffusivity
t_end = 0.001; % final time
% Grid parameters
[nx,ny] = size(xx);
dx = xx(2,1)-xx(1,1);
dy = yy(1,2)-yy(1,1);
% determine time step from stability requirenment
safety_fac = 1.1; % safety factor
dt_max = 0.25/(1/dx^2+1/dy^2)/alpha; % maximum time step
nt = ceil(t_end/(dt_max/safety_fac)); % number of time steps
dt = t_end/nt; % time step that gets us right to
t_end
disp(['Stability requirenment: ' num2str(dt*alpha*(1/dx^2+1/dy^2)) ' < 0.25'])

figure
for ti = 1:nt
%%%%%%%%%%%%%%%%%%
% EXPLICIT EULER % (1st-order forward difference in time)
%%%%%%%%%%%%%%%%%%
d2Tdx2 = d2dx2(T,dx,'periodic');
d2Tdy2 = d2dy2(T,dy,'periodic');
T = T + alpha*dt*d2Tdx2 + alpha*dt*d2Tdy2;
% Plot temperature at time t
subplot(1,2,1)
pcolor(xx,yy,T); shading interp; axis equal tight;
colorbar;
title(['T @ t=' num2str(dt*ti)])
xlabel('x');
ylabel('y');
caxis([0 1])
subplot(1,2,2)
plot(yy(1,:),T(80,:)), ylim([min(T(80,:)),max(T(80,:))])
drawnow
end
%%%%%%%%%%%%%%%%%%
% FUNCTIONS %
%%%%%%%%%%%%%%%%%%
% Second-order second central difference function in x
function d2fdx2 = d2dx2(f,dx,bc)
% set default value for 'bc'
if nargin<3, bc = 'one-sided'; end
% determine field size
[nx,ny] = size(f);
% allocate return field
d2fdx2 = zeros(nx,ny);
dx2 = dx^2;
% central difference
for i=2:nx-1
for j=1:ny
d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(i-1,j))/dx2;
end
end
switch bc
case 'periodic'
i = 1;
for j=1:ny
d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(nx,j))/dx2;
end
i = nx;
for j=1:ny
d2fdx2(i,j) = (f(1,j)-2*f(i,j)+f(i-1,j))/dx2;
end
otherwise
% forward difference for first point
i = 1;
for j=1:ny
d2fdx2(i,j) = (2*f(i,j)-5*f(i+1,j)+4*f(i+2,j)-f(i+3,j))/dx2;
end
% backward difference for last point
i = nx;
for j=1:ny
d2fdx2(i,j) = (2*f(i,j)-5*f(i-1,j)+4*f(i-2,j)-f(i-3,j))/dx2;
end
end
end
% Second-order second central difference function in y
function d2fdy2 = d2dy2(f,dy,bc)
% set default value for 'bc'
if nargin<3, bc = 'one-sided'; end
% determine field size
[nx,ny] = size(f);
% allocate return field
d2fdy2 = zeros(nx,ny);
dy2 = dy^2;
% central difference
for i=1:nx
for j=2:ny-1
d2fdy2(i,j) = (f(i,j+1)-2*f(i,j)+f(i,j-1))/dy2;
end
end
switch bc
case 'periodic'
j = 1;
for i=1:nx
d2fdy2(i,j) = (f(i,j+1)-2*f(i,j)+f(i,ny))/dy2;
end
j = ny;
for i=1:nx
d2fdy2(i,j) = (f(i,1)-2*f(i,j)+f(i,j-1))/dy2;
end
otherwise
% forward difference for first point (bottom)
j = 1;
for i=1:nx
d2fdy2(i,j) = (2*f(i,j)-5*f(i,j+1)+4*f(i,j+2)-f(i,j+3))/dy2;
end
% backward difference for last point (top)
j = ny;
for i=1:nx
d2fdy2(i,j) = (2*f(i,j)-5*f(i,j-1)+4*f(i,j-2)-f(i,j-3))/dy2;
end
end
end


