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
dfdy_analytical     = -cos(xx).*sin(yy);
d2fdy2_analytical   = -cos(xx).*cos(yy);

%% First derivative and error
figure
subplot(2,5,[1 6])
pcolor(xx,yy,f)
caxis([-1 1])
axis equal tight
title('f')

subplot(2,5,[2 7])
pcolor(xx,yy,dfdy_analytical)
caxis([-1 1])
axis equal tight
title('df/dy (analytical)')

subplot(2,5,3)
dfdy = ddy_fwd(f,dy);
pcolor(xx,yy,dfdy)
caxis([-1 1])
axis equal tight
title('df/dy (forward difference)')

subplot(2,5,3+5)
pcolor(xx,yy,abs(dfdy-dfdy_analytical))
caxis([-0.1 0.1])
axis equal tight
title('Error - df/dy (forward difference)')

subplot(2,5,4)
dfdy = ddy_bwd(f,dy);
pcolor(xx,yy,dfdy)
caxis([-1 1])
axis equal tight
title('df/dy (backward difference)')

subplot(2,5,4+5)
pcolor(xx,yy,abs(dfdy-dfdy_analytical))
caxis([-0.1 0.1])
axis equal tight
title('Error - df/dy (backward difference)')


subplot(2,5,5)
dfdy = ddy_central(f,dy);
pcolor(xx,yy,dfdy)
caxis([-1 1])
axis equal tight
title('df/dy (central difference)')

subplot(2,5,5+5)
pcolor(xx,yy,abs(dfdy-dfdy_analytical))
caxis([-0.1 0.1])
axis equal tight
title('Error - df/dy (central difference)')







