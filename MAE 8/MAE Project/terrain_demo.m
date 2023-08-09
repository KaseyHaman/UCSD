

load('terrain.mat');

%% Estimate the terrain elevation h at Xo = 70.5 m and Yo = 20.5 m
Xo = 70.5;
Yo = 20.5;
h = interp2(x_terrain,y_terrain,h_terrain,Xo,Yo);


%% Plot terrain surface
figure(1); hold on;
plot3(X,Y,Z,'-r', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Terrain Demo');
legend('target location on terrain','terrain surface')
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
    'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);