%% Test scripts for finite difference functions

clear all
close all
clc

% make sure to copy your *.m-files into the 'functions' directory  
addpath('functions')

% 1-D coordinates 
nx          = 50;
ny          = 30;
x           = linspace(0,10,nx);
y           = linspace(0,7, ny);
dx          = x(2)-x(1);
dy          = y(2)-y(1);

% 2-D coordinates
[xx,yy]     = ndgrid(x,y);

% surrogate data
f                   =  cos(xx).*cos(yy);
dfdx_analytical     = -sin(xx).*cos(yy);
d2fdx2_analytical   = -cos(xx).*cos(yy);

%% Second derivative and error
figure
subplot(1,4,1)
pcolor(xx,yy,f)
caxis([-1 1])
axis equal tight
title('f (analytical)')

subplot(1,4,2)
pcolor(xx,yy,d2fdx2_analytical)
caxis([-1 1])
axis equal tight
title('d^2f/dx^2 (analytical)')

subplot(1,4,3)
d2fdx2 = d2dx2(f,dx);
pcolor(xx,yy,d2fdx2)
caxis([-1 1])
axis equal tight
title('d^2f/dx^2 (central second difference)')

subplot(1,4,4)
pcolor(xx,yy,abs(d2fdx2-d2fdx2_analytical))
caxis([-0.1 0.1])
axis equal tight
title('Error - d^2f/dx^2 (central second difference)')





