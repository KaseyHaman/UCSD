%% MAE 171 Fluids Lab
clc; clear all; close all;

% MATLAB Script: Plot Uniform Flow, Rankine Full Body, and Half Body with Obstacle

% Parameters
U = 1.0;                  % Freestream velocity (m/s)
mu = 1.0;                 % Doublet strength (m^2/s)
Q_s = 3.0;                % Source strength (m^3/s)
Q_sink = 3.0;             % Sink strength (m^3/s)
circle_diameter = 5* 4.5 / 39.37; % Convert 4.5 inches to meters
x_min = -2;               % x-axis lower limit (m)
x_max = 4;                % x-axis upper limit (m)
y_min = -2;               % y-axis lower limit (m)
y_max = 2;                % y-axis upper limit (m)
x_source = 0;             % x-coordinate of the source (m)
y_source = 0;             % y-coordinate of the source (m)
x_sink = 1;               % x-coordinate of the sink (m)
y_sink = 0;               % y-coordinate of the sink (m)
x_cen = 0; y_cen=0;       % Center values

% Create Grid
x = linspace(x_min, x_max, 500);  % High grid resolution for smooth plot
y = linspace(y_min, y_max, 500);
[X, Y] = meshgrid(x, y);          % Create mesh grid for flow field

% Compute Stream Function for Uniform Flow
psi_uniform = U * Y;  % Uniform flow contribution

% Compute Stream Function for Source
psi_source = (Q_s / (2 * pi)) * atan2(Y - y_source, X - x_source); % Source contribution

% Compute Stream Function for Sink
psi_sink = (Q_sink / (2 * pi)) * atan2(Y - y_sink, X - x_sink);   % Sink contribution

% Compute Total Stream Function for Rankine Full Body
psi_rankine = psi_uniform + psi_source - psi_sink; % Rankine full body

%Compite Rankine half body
psi_half_body = psi_uniform + psi_source; % Rankine half-body (no sink)

%Now we need to compute flow around a cylinder

%Define the cylinder
theta = linspace(0, 2*pi, 100);
x_circle = circle_diameter / 2 * cos(theta);  % Circle boundary
y_circle = circle_diameter / 2 * sin(theta);  % Circle boundary

% Define cylinder radius
R = circle_diameter / 2;

% Update doublet strength
mu = 2 * pi * U * R^2;

% Compute Stream Function for Doublet Flow
psi_doublet = -mu * Y ./ (2 * pi * (X.^2 + Y.^2));

% Total Stream Function for Flow Around Cylinder
psi_cylinder = psi_uniform + psi_doublet;

% Plot Streamlines for Flow Around a Circular Cylinder
figure;
contour(X, Y, psi_cylinder, 100);    % Plot total flow streamlines
hold on;
plot(x_circle, y_circle, 'r-', 'LineWidth', 2); % Plot the circle boundary (obstacle)
title('Flow Around Circular Cylinder');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
grid on;
hold off;

% Plot Streamlines for Uniform Flow
figure;
contour(X, Y, psi_uniform, 30, 'b');    % Plot uniform flow streamlines
title('Uniform Flow Analytical Solution');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
grid on;
hold off;

% Plot Streamlines for Rankine Full Body Flow
figure;
contour(X, Y, psi_rankine, 30, 'g');    % Plot Rankine full body streamlines
title('Rankine Full Body Flow Analytical Simulation');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
grid on;
hold off;

% Plot Streamlines for Rankine Half Body Flow
figure;
contour(X, Y, psi_half_body, 30, 'm');    % Plot Rankine half body streamlines
title('Rankine Half Body Analytical Simulation');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
grid on;
hold off;




% 
% % Parameters
% U =0.00008832;                  % Freestream velocity (m/s)
% mu = 1.0;                 % Doublet strength (m^2/s)
% Q_s = 0.000003469;                % Source strength (m^3/s)
% Q_sink = 0.000003469*4/3;             % Sink strength (m^3/s)
% circle_diameter = 4.5 / 39.37; % Convert 4.5 inches to meters
% Qmodel_diameter = 1 / 39.37;  % Conver 1.5 inches to meters
% x_min = -0.455;               % x-axis lower limit (m)
% x_max = 0.455;                % x-axis upper limit (m)
% y_min = -0.2925;               % y-axis lower limit (m)
% y_max = 0.2925;                % y-axis upper limit (m)
% x_source = x_min + (8.5 / 39.37);             % x-coordinate of the source (m)
% y_source = 0;             % y-coordinate of the source (m)
% x_sink = x_source + (8 / 39.37);               % x-coordinate of the sink (m)
% y_sink = 0;               % y-coordinate of the sink (m)
% x_cen = 0; y_cen=0;       % Center values