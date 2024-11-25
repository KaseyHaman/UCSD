clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%% Problem 1

%Part A
sigma = 1.22;
x = linspace(0, 10, 10000); 
f_x = (x / sigma^2) .* exp(-(x.^2) / (2 * sigma^2)); % PDF Equation

figure(1);
plot(x, f_x, '-b', 'LineWidth', 1.5);
xlabel('Wind Speed');
ylabel('PDF Output');
title('Rayleigh Distribution of Wind Turbines');
hold on;

%Parts B-D Done on paper

%Part E
N = 1000000;
r = rand(N, 1);
x_r = sigma*sqrt(-2*log(1-r));% Inputting r (= percentile in distribution)
% to find which x value corresponds to that percentile in the PDF.
% GUESS: The area of the CDF is already 1, so adding up all the x values
% and dividing by the number of x values gives your most expected x value.

figure(2);
histogram(x_r, 'Normalization', 'pdf')
hold on;
plot(x, f_x, 'c', 'LineWidth', 2);
title('Experimental Rayleigh Distribution');
xlabel('Wind Speed');
ylabel('PDF Output');

%Part F
mean(x_r); %The value is approximately equal to the expected value of x that was calculated.

%Part G
x_raylrnd = raylrnd(sigma, N, 1);
figure(3);
histogram(x_raylrnd, 'Normalization', 'pdf');
hold on;
plot(x, f_x, 'c', 'LineWidth', 2);
xlabel('Wind Speed');
ylabel('PDF Output');
title('MATLAB Rayleigh Distribution');



%% Problem 2
%Part A
clear; clc;
x = linspace(-1000, 1000, 10000); 
% Variables
m = 6.63*(10^(-26));
T = 273;
k = (1.38*10^-23);
% PDF Equation
f_x = (sqrt(m/(2*pi*k*T))*exp(-(m*x.^2)/(2*k*T))); 
figure(4);
plot(x, f_x, '-b', 'LineWidth', 1.5);
xlabel('Velocity Projection');
ylabel('PDF Output');
title('Argon Velocity Distribution');
hold on;

%Part B Calculations done on paper

v = linspace(0, 2000, 1000);
g_v = (m / (2 * pi * k * T))^(3/2) .* 4 * pi * v.^2 .* exp(-m * v.^2 / (2 * k * T));
figure(5);
plot(v, g_v);
xlabel('Speed');
ylabel('g(v)');
title('Distribution of Absolute Velocities for Argon');

%Part C
N = 1000000;
vx = randn(N, 1) * sqrt(k * T / m);
vy = randn(N, 1) * sqrt(k * T / m);
vz = randn(N, 1) * sqrt(k * T / m);
v = sqrt(vx.^2 + vy.^2 + vz.^2); %Note sigma = sqrt(kt/m) by gaussian
%equation. So we use this to determine our v values.
figure(6);
histogram(v, 100);
title('Simulated Absolute Velocities');
xlabel('Speed');
ylabel('PDF Output');

figure(7);
histogram(vx, 100);
title('Simulated Vx Velocity');
xlabel('Speed');
ylabel('PDF Output');


%% Problem 3
%Part A
N = 1000000;
F = 300 + 10 * (2 * rand(N, 1) - 1);
D = 0.04 + 0.0001 * (2 * rand(N, 1) - 1);
d = 0.004 + 0.0004 * (2 * rand(N, 1) - 1);
figure(8)
histogram(F, 100)
title('Randomly Distributed F');
figure(9)
histogram(D, 100)
title('Randomly Distributed D');
figure(10)
histogram(d, 100)
title('Randomly Distributed d');

%Part B
C = D ./ d;
K = (4 * C - 1) ./ (4 * C - 4) + 0.615 ./ C;
tau_max = K .* (8 * F .* D) ./ (pi * d.^3);
figure(11)
histogram(tau_max, 100)
title('Maximum Shear Stress From Random Sampling');

%Part C
clear F; clear D; clear d; clear C; clear K; clear tau_max;
F = normrnd(300, 10/3, [N, 1]);
D = normrnd(0.04, 0.0001/3, [N, 1]);
d = normrnd(0.004, 0.0004/3, [N, 1]);
figure(12)
histogram(F, 100)
title('Normally Distributed F');
figure(13)
histogram(D, 100)
title('Normally Distributed D');
figure(14)
histogram(d, 100)
title('Normally Distributed d');

% Repeat stress calculations as above
C = D ./ d;
K = (4 * C - 1) ./ (4 * C - 4) + 0.615 ./ C;
tau_max = K .* (8 * F .* D) ./ (pi * d.^3);
figure(15)
histogram(tau_max, 100)
title('Maximum Shear Stress Deviation');


