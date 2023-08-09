clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';
hw_num = 'project';
%%x, y, z = position (meters)
%%t = time
%%U, V, W = velocity components (m/s in x, y and z direction respectively)
global x_terrain y_terrain h_terrain
load('terrain.mat');
%%Task 1
exp_num = 1;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eom1 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eom1mass = m;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(1); subplot(2, 3, 1);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z, '-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast')
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 2;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eom2 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eom2mass = m;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(1); subplot(2, 3, 2);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 3;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eom3 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eom3mass = m;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(1); subplot(2, 3, 3);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 4;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eom4 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eom4mass = m;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(1); subplot(2, 3, 4);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 5;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
Eom5 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eom5mass = m;
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(1); subplot(2, 3, 5);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 6;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eom6 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eom6mass = m;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(1); subplot(2, 3, 6);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 7;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eos7 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eos7spring = k;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(2); subplot(2, 3, 1);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z, '-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast')
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 8;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eos8 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eos8spring = k;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
figure(2); subplot(2, 3, 2);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 9;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eos9 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eos9spring = k;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(2); subplot(2, 3, 3);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 10;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eos10 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eos10spring = k;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(2); subplot(2, 3, 4);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 11;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eos11 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eos11spring = k;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(2); subplot(2, 3, 5);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 12;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eos12 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eos12spring = k;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(2); subplot(2, 3, 6);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 13;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eoc13 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eoc13cord = l;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(3); subplot(2, 3, 1);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z, '-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast')
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 14;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eoc14 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eoc14cord = l;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(3); subplot(2, 3, 2);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 15;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eoc15 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eoc15cord = l;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(3); subplot(2, 3, 3);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 16;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eoc16 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eoc16cord = l;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(3); subplot(2, 3, 4);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 17;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eoc17 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eoc17cord = l;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(3); subplot(2, 3, 5);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);
exp_num = 18;
read_input('bungee_data.txt', exp_num);
global m k l Xo Yo Zo Uo Vo Wo;
[T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo);
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = [];
Eoc18 = max(sqrt(((X(:)).^2 + (Y(:)).^2 + (Z(:)).^2)));
Eoc18cord = l;
expsafety(exp_num) = safety;
global Vmag
vmag(exp_num).values = Vmag;
mass(exp_num).value = m;
x(exp_num).values = X;
y(exp_num).values = Y;
z(exp_num).values = Z;
if safety == 1
    safetystr = [];
else
    safetystr = 'danger';
end
strtitle = sprintf('Exp. #%1i %s\n', exp_num, safetystr);
figure(3); subplot(2, 3, 6);
plot3(X(end), Y(end), Z(end), 'ro', 'Markerfacecolor', 'r', 'Markersize', 12)
hold on;
plot3(X,Y,Z,'-k', 'linewidth', 2);
surf(x_terrain, y_terrain, h_terrain);
colormap('jet'); shading interp;
grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(strtitle);
% legend('Bungee Jump Path','Terrain Surface', 'Location', 'northeast', 'fontsize', 4)
view(3); axis([-100 100 -150 100 -250 0]);
set(gca, 'LineWidth', 2, 'FontSize', 10, ...
  'Xtick',-150:50:150, 'Ytick', -150:50:100, 'Ztick',-200:50:0);

for n = 1:18
    if expsafety(n) == 1
        str(n).s = sprintf('Exp #%i', n);
    else
        str(n).s = sprintf('Exp #%i danger', n);
    end
end

xdata = [Eom1mass; Eom2mass; Eom3mass; Eom4mass; Eom5mass; Eom6mass];
ydata = [Eom1; Eom2; Eom3; Eom4; Eom5; Eom6];
figure(4); subplot(1, 3, 1); hold on;
plot (Eom1mass, Eom1, 'ko', 'Markersize', 12, 'Markerfacecolor', 'k')
plot (Eom2mass, Eom2, 'ro', 'Markersize', 12, 'Markerfacecolor', 'r')
plot (Eom3mass, Eom3, 'bo', 'Markersize', 12, 'Markerfacecolor', 'b')
plot (Eom4mass, Eom4, 'go', 'Markersize', 12, 'Markerfacecolor', 'g')
plot (Eom5mass, Eom5, 'co', 'Markersize', 12, 'Markerfacecolor', 'c')
plot (Eom6mass, Eom6, 'mo', 'Markersize', 12, 'Markerfacecolor', 'm')
plot(xdata, ydata, '-k')
grid on;
xlabel('Mass (kg)');
ylabel('Max Distance (m)');
title('Effect of Mass');
legend(str(1).s, str(2).s, str(3).s, str(4).s, str(5).s, str(6).s, 'location', 'best')
xdata = [Eos7spring; Eos8spring; Eos9spring; Eos10spring; Eos11spring; Eos12spring];
ydata = [Eos7; Eos8; Eos9; Eos10; Eos11; Eos12];
axis('tight');
figure(4); subplot(1, 3, 2); hold on;
plot (Eos7spring, Eos7, 'ko', 'Markersize', 12, 'Markerfacecolor', 'k')
plot (Eos8spring, Eos8, 'ro', 'Markersize', 12, 'Markerfacecolor', 'r')
plot (Eos9spring, Eos9, 'bo', 'Markersize', 12, 'Markerfacecolor', 'b')
plot (Eos10spring, Eos10, 'go', 'Markersize', 12, 'Markerfacecolor', 'g')
plot (Eos11spring, Eos11, 'co', 'Markersize', 12, 'Markerfacecolor', 'c')
plot (Eos12spring, Eos12, 'mo', 'Markersize', 12, 'Markerfacecolor', 'm')
plot(xdata, ydata, '-k')
grid on;
xlabel('Spring Modulus (N/m)');
ylabel('Max Distance (m)');
title('Effect of Spring Modulus');
legend(str(7).s, str(8).s, str(9).s, str(10).s', str(11).s, str(12).s, 'location', 'best')
xdata = [Eoc13cord; Eoc14cord; Eoc15cord; Eoc16cord; Eoc17cord; Eoc18cord];
ydata = [Eoc13; Eoc14; Eoc15; Eoc16; Eoc17; Eoc18];
axis('tight');
figure(4); subplot(1, 3, 3); hold on;
plot (Eoc13cord, Eoc13, 'ko', 'Markersize', 12, 'Markerfacecolor', 'k')
plot (Eoc14cord, Eoc14, 'ro', 'Markersize', 12, 'Markerfacecolor', 'r')
plot (Eoc15cord, Eoc15, 'bo', 'Markersize', 12, 'Markerfacecolor', 'b')
plot (Eoc16cord, Eoc16, 'go', 'Markersize', 12, 'Markerfacecolor', 'g')
plot (Eoc17cord, Eoc17, 'co', 'Markersize', 12, 'Markerfacecolor', 'c')
plot (Eoc18cord, Eoc18, 'mo', 'Markersize', 12, 'Markerfacecolor', 'm')
plot(xdata, ydata, '-k')
grid on;
xlabel('Cord Length (m)');
ylabel('Max Distance (m)');
title('Effect of Cord Length');
legend(str(13).s, str(14).s, str(15).s, str(16).s, str(17).s, str(18).s, 'location', 'best')
axis('tight');
%%Task 2
% mass(exp_num).value = m;
% x(exp_num).values = X;
% y(exp_num).values = Y;
% z(exp_num).values = Z;
dt = 0.02;
for n = 1:18
  exp_res(n).number = n;
  exp_res(n).max_speed = max(vmag(n).values);
  exp_res(n).max_acceleration = max(diff((vmag(n).values)./dt));
  KE(n).values = 0.5.*mass(n).value.*(vmag(n).values).^2;
  exp_res(n).integrated_KE = dt * (0.5 * (KE(n).values(end)+ KE(n).values(1)) + sum(KE(n).values(2:end-1)));
  exp_res(n).travel_distance = sum(sqrt((diff(x(n).values)).^2 + (diff(y(n).values)).^2 + (diff(z(n).values)).^2));
end
%%Task 3
fid = fopen('report.txt', 'wt');
fprintf(fid, 'Kasey Haman\nA16978114\n');
fprintf(fid, 'exp number, max speed (m/s), max acc (m/sˆ2), int KE (J s), travel dist (m)\n');
for n = 1:18
fprintf(fid, '%10d %15.7e %15.7e %15.7e %15.7e\n', exp_res(n).number, exp_res(n).max_speed,...
  exp_res(n).max_acceleration, exp_res(n).integrated_KE, exp_res(n).travel_distance);
end
fclose(fid);
safety = expsafety(:)';
p1a ='See figure 1';
p1b ='See figure 2';
p1c ='See figure 3';
p1d ='See figure 4';
p2a = exp_res(1);
p2b = [exp_res.max_speed];
p2c = [exp_res.max_acceleration];
p2d = [exp_res.integrated_KE];
p2e = [exp_res.travel_distance];
p3 = evalc('type report.txt');

