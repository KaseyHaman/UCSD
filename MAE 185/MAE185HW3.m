%% 3.1
clc;clear all;
%Load data
load('gotritons','T','xx','yy')

%Define Step Sizes
dx = xx(2, 1)-xx(1, 1);
dy = yy(1, 2)-yy(1, 1);

%Define variables
alpha = 2;
SF = 2;
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
    drawnow;

    T = T + dt * (-ddx_central_periodic(T, dx) - ddy_central_periodic(T, dy)); %Set new T

    t_central = t_central + dt; %increase time

   
end

