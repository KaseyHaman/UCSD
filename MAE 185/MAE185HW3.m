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