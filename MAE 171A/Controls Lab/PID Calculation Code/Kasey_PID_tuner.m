clc;
clear;
close all;
%% 
% Define the transfer function G(s)
num = [6.04e-09, 1.088e-08, 1.434e-06];
den = [2.05e-15, 4.409e-13, 4.2e-12, 3.369e-10, 1.002e-09, 3.253e-08];
G = tf(num, den);

% Step Response of Open-Loop System
figure;
step(G);
title('Open-Loop Step Response');
grid on;

% Force a full PID controller (use 'pidstd' structure)
C = pidtune(G, 'PID', 50); % 50 rad/s is an initial crossover frequency guess

% Display PID Gains
disp('Initial PID Gains:');
disp(C);

% Closed-Loop System with PID Controller
T = feedback(C*G, 1);

% Step Response of the Closed-Loop System
figure;
step(T);
title('Closed-Loop Step Response with PID');
grid on;

% Extract Performance Metrics
info = stepinfo(T);

% Display Performance Results
disp('Closed-Loop Performance:');
fprintf('Overshoot: %.2f%%\n', info.Overshoot);
fprintf('Settling Time: %.4f seconds\n', info.SettlingTime);


%%
% Manual Fine-Tuning for Better Performance
Kp = C.Kp;
Ki = C.Ki;
Kd = C.Kd;

% Adjust gains for better transient response
Kp_new = Kp * 1.2;  % Increase proportional gain for faster response
Ki_new = Ki * 1.1;  % Slightly increase integral gain for better steady-state
Kd_new = Kd * 0.7;  % Reduce derivative to avoid excessive damping

C_tuned = pid(Kp_new, Ki_new, Kd_new);

% New Closed-Loop System with Tuned PID
T_tuned = feedback(C_tuned*G, 1);

% Step Response of the Tuned System
figure;
step(T_tuned);
title('Tuned Closed-Loop Step Response');
grid on;

% Final Performance Metrics
info_tuned = stepinfo(T_tuned);
disp('Tuned Closed-Loop Performance:');
fprintf('Overshoot: %.2f%%\n', info_tuned.Overshoot);
fprintf('Settling Time: %.4f seconds\n', info_tuned.SettlingTime);

% Bode Plot Analysis to Check Stability Margins
figure;
margin(C_tuned * G);
title('Bode Plot of Tuned PID Controller');
grid on;

%Check ss error
ess = 1 - dcgain(T_tuned);
disp(['Steady-State Error: ', num2str(ess)]);

