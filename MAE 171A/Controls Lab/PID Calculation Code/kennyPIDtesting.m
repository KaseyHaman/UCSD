clear
close all
clc

% Given Transfer Function
num = [6.04e-09 1.088e-08 1.434e-06];
den = [2.05e-15 4.409e-13 4.2e-12 3.369e-10 1.002e-09 3.253e-08];

% Define the system transfer function
G = tf(num, den);

% Define PID controller (adjusting the parameters as needed)
% Kp = 0.0663;  % Proportional gain
% Ki = 0;  % Integral gain
% Kd = 0;  % Derivative gain

% [Kp, Ki, Kd] = ZieglerNichols(0.0663, 0.17, 'ClassicPID');

Kp = 0.2;
Ki = 0.468;
Kd = 0.02;

% Create the PID controller
C = pid(Kp, Ki, Kd);

% Closed-loop system with PID controller
sys_closed_loop = feedback(C*G, 1);

% scaling for steady state value
scale = 44;

% Plot both the open-loop and closed-loop step responses
figure;
subplot(2,1,1);
step(G);
title('Original Transfer Function Step Response');
grid on;
yline(1.25*scale, '--b', '+25% (Overshoot Limit)', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'top')
yline(0.98*scale, '--r', '-2% (Steady State Min)', 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment', 'bottom');
yline(1.02*scale, '--g', '+2% (Steady State Max)', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'top');
xlim([0 30]);  % Set x-axis to 10 seconds

subplot(2,1,2);
step(scale*sys_closed_loop);
title('PID-Controlled System Step Response');
grid on;
hold on
yline(1.25*scale, '--b', '+25% (Overshoot Limit)', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'top')
yline(0.98*scale, '--r', '-2% (Steady State Min)', 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment', 'bottom');
yline(1.02*scale, '--g', '+2% (Steady State Max)', 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'top');
xlim([0 30]);  % Set x-axis to 10 seconds

