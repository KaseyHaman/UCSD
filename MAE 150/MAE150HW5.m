clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%% Problem 3
%Load/Visualize Data
TP = readmatrix('cam_contour.dat');
x = TP(:, 1);
y = TP(:, 2);
% 
% xt = [x(380:end); x(1:30)];
% yt = [y(380:end); y(1:30)];
figure(1);
plot(x,y, '-b', 'LineWidth', 1);
grid on;
hold on;
axis on;
xlabel('x')
ylabel('y')

%ATTEMPT 1
%Gradient function is more accurate than diff function
% dx = gradient(x);
% dy = gradient(y);
% d2x = gradient(dx);
% d2y = gradient(dy);
% 
% %Formula for p
% p = ((dx).^2+(dy).^2).^(3/2)./abs(dx.*d2y-d2x.*dy);

% %ATTEMPT 2
r_data = sqrt(x.^2 + y.^2);
theta_data = atan2(y, x);
d_r_data = gradient(r_data);
d_r_data2 = gradient(d_r_data);

p = ((r_data.^2+d_r_data.^2).^(3/2)./(r_data.^2+2*(d_r_data.^2)-r_data.*d_r_data2));


%The radius of the roller must be less than the min p value
print('The roller must have a radius less than')
min(p)


%% Problem 4
%% Part A
%Done on paper

%% Part B
L = 5; %mm 
riset = 0.02; %s
dwellt = 0.1; %s
fallt = 0.02; %s
omega = 2*pi / (riset+dwellt+fallt);

trise = linspace(0, riset, 100);
tdwell = linspace(riset, riset+dwellt, 100);
tfall = linspace(riset+dwellt, riset+dwellt+fallt, 100);

%Solve rise time case
beta_rise = omega*riset; %Radian amount that rise time goes through
thetarise = omega .* trise;

yrise = L * (thetarise / beta_rise - (1 / (2 * pi)) * sin(2 * pi * thetarise / beta_rise));
vrise = L / beta_rise * (1 - cos(2 * pi * thetarise / beta_rise)) * omega; 
arise = L * (2 * pi / beta_rise^2) * sin(2 * pi * thetarise / beta_rise) * omega^2; 

%Solve dwell time case
beta_dwell = omega*dwellt; %Radian amount that rise time goes through
thetadwell = omega .* tdwell;
ydwell = L*ones(1, size(thetadwell, 2));
vdwell = zeros(1, size(thetadwell, 2));
adwell = zeros(1, size(thetadwell, 2));

%Solve fall time case
beta_fall = omega*fallt; %Radian amount that rise time goes through
thetafall = omega .* tfall;

yfall = L - L * (thetafall - thetafall(1)) / beta_fall + (L / (2 * pi)) * sin(2 * pi * (thetafall - thetafall(1)) / beta_fall);
vfall = -L / beta_fall * (1 - cos(2 * pi * thetafall / beta_fall)) * omega; 
afall = -L * (2 * pi / beta_fall^2) * sin(2 * pi * thetafall / beta_fall) * omega^2; 

%Consolidate data
t = [trise, tdwell, tfall];
x = [yrise, ydwell, yfall];
v = [vrise, vdwell, vfall];
a = [arise, adwell, afall];
theta = [thetarise, thetadwell, thetafall];


%Plot data
figure (2);
subplot(3, 1, 1);
plot (t, x, '-b')
hold on;
title('Displacement')
ylabel('y');
xlabel('t');

subplot(3, 1, 2);
plot(t, v,'-k')
title('Velocity')
ylabel('v');
xlabel('t');

subplot(3, 1, 3);
plot(t, a,'-r')
title('Acceleration')
ylabel('a');
xlabel('t');


%% Part C 
%Given values
L = 5; %mm 
r = 2; %mm
phi = pi/6; %rad
riset = 0.02; %s
dwellt = 0.1; %s
fallt = 0.02; %s
omega = 2*pi / (riset+dwellt+fallt); 
beta = omega * riset;
gamma = omega * fallt;

%Max velocity
Vmax = max(vrise);

%Find Rb
Rprime = r + L; % Prime circle radius simplification
for i = 1:length(theta)
    dy_dtheta = L / beta_rise * (1 - cos(2 * pi * theta(i) / beta_rise)); % dy/dtheta from rise profile
    current_phi = atan((dy_dtheta / omega) / (Rprime + y(i))); % Current pressure angle
    if current_phi <= phi
        Rb = Rprime - y(i); % Update Rb
    end
end

% Calculate y(theta) for rise, dwell, and fall
y_rise = L * (thetarise / (omega * riset) - (1 / (2 * pi)) * sin(2 * pi * thetarise / (omega * riset)));
y_dwell = L * ones(size(thetadwell));
y_fall = L * (1 - (thetafall - thetadwell(end)) / (omega * fallt) + (1 / (2 * pi)) * sin(2 * pi * (thetafall - thetadwell(end)) / (omega * fallt)));

%Define theta and y vectors
y = x;

% Cam and pitch contours
cam_contour_x = (Rb + y) .* cos(theta);
cam_contour_y = (Rb + y) .* sin(theta);
pitch_curve_x = (Rb + y + r) .* cos(theta);
pitch_curve_y = (Rb + y + r) .* sin(theta);

%% Part D 
% Pitch point for velocity max
[~, pitch_index] = max(y);
pitch_point_x = pitch_curve_x(pitch_index);
pitch_point_y = pitch_curve_y(pitch_index);
pitch_point_radius = sqrt(pitch_point_x^2+pitch_point_y^2);

% Plot circles and curves
figure(3);
theta_circle = linspace(0, 2 * pi, 200);
x_base = Rb * cos(theta_circle);
y_base = Rb * sin(theta_circle);
x_prime = (Rb + r) * cos(theta_circle);
y_prime = (Rb + r) * sin(theta_circle);

hold on;
plot(x_base, y_base, 'cyan', 'LineWidth', 2);
plot(x_prime, y_prime, 'green', 'LineWidth', 2);
plot(pitch_curve_x, pitch_curve_y, 'blue', 'LineWidth', 2);
plot(cam_contour_x, cam_contour_y, 'black', 'LineWidth', 2);
plot(pitch_point_x, pitch_point_y, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'red');
plot(pitch_point_radius * cos(theta_circle), pitch_point_radius * sin(theta_circle), 'r--');
title('Cam Profile with Pitch Point and Circle');
legend('Base Circle', 'Prime Circle', 'Pitch Curve', 'Cam Contour', 'Pitch Point', 'Pitch Circle');
axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
grid on;
hold off;

% %% Part E
figure(4);
axis equal;
grid on;
hold on;
xlabel('X (mm)');
ylabel('Y (mm)');
title('Cam-Follower Mechanism Animation');
plot(x_base, y_base, 'cyan', 'LineWidth', 2); % Base circle
  pause(0.1);
for k = 1:length(theta)
    cam_x = (Rb + y(k)) * cos(theta(k));
    cam_y = (Rb + y(k)) * sin(theta(k));
    follower_y = y(k) + r;

    cla;
    plot(x_base, y_base, 'cyan', 'LineWidth', 2);
    plot(cam_contour_x, cam_contour_y, 'black', 'LineWidth', 1);
    plot(pitch_curve_x, pitch_curve_y, 'blue', 'LineWidth', 1);

    plot([0 cam_x], [0 cam_y], 'k', 'LineWidth', 2);
    plot(cam_x, cam_y, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'black');
    plot([0 0], [0 follower_y], 'r', 'LineWidth', 3);
    plot(0, follower_y, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');

    pause(0.05);
end

hold off;