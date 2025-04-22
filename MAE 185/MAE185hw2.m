clc;clear all;


%Load data
load('cylinder_Re100.mat')

%Size u
[a b c] = size(u);
dx = x(2)-x(1);
dy = y(1, 2)-y(1, 1);


%Run simulation
for t = 1:a
figure(1);
    ustep = squeeze(u(t, :, :));
    dudy = ddy_central(ustep, dy);

    vstep = squeeze(v(t, :, :));
    dvdx = ddx_central(vstep, dx);

    vorticity = dvdx - dudy;

    pcolor(x, y, vorticity);
    shading interp;
    axis equal tight;
    title('Vorticity');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    drawnow;

end