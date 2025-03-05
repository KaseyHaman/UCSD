clear all;
close all;
clc;

%% Problem 2
% Given tensile strength data (MPa)
sigma = [321, 389, 411, 423, 438, 454, 475, 489, 497, 501];

% Sort the data in ascending order
sigma = sort(sigma);

% Number of data points
N = length(sigma);

% Weibull probability estimator
for i = 1:N
    num = N - i;
P(i) = (num+1) ./ (N + 1);
end
% Weibull plot transformation
ln_sigma = log(sigma);
ln_ln_invP = log(log(1./P));

% Perform linear regression to find slope (Weibull modulus)
coeffs = polyfit(ln_sigma, ln_ln_invP, 1);
m = coeffs(1); % Weibull modulus
C = exp(-coeffs(2)/m); % Scale parameter

% Plot Weibull distribution
figure;
plot(ln_sigma, ln_ln_invP, 'o', 'MarkerFaceColor', 'b'); % Data points
hold on;
plot(ln_sigma, polyval(coeffs, ln_sigma), 'r-', 'LineWidth', 2); % Fit line
grid on;
xlabel('ln(sigma)');
ylabel('ln(ln(1/P(V)))');
title('Weibull Plot Problem 2');
legend('Data Points', 'Linear Fit', 'Location', 'best');

% Scaling for 60 mm specimen using length effect
L_ref = 20; % Reference length (mm)
L_new = 60; % New length (mm)
sigma_50 = C * (log(2))^(1/m) * (L_ref / L_new)^(1/m); % 50% survival probability

% Display results
fprintf('Weibull Modulus (m): %.2f\n', m);
fprintf('Scale Parameter (C): %.2f MPa\n', C);
fprintf('Tensile Strength at 50%% survival for 60mm specimen: %.2f MPa\n', sigma_50);


%% Problem 3


% Given failure loads (N)
P_loads = [1040, 1092, 1120, 1210, 1320, 1381, 1419, 1470, 1490, 1540];

% Beam dimensions (mm)
b = 10; % Width
h = 5;  % Height

% Span length (mm)
L = 50;

% Calculate flexural strength using sigma = 3PL / (2bh^2)
sigma = (3 .* P_loads .* L) ./ (2 * b * h^2); % MPa

% Sort the data in ascending order
sigma = sort(sigma);

% Number of data points
N = length(sigma);

% Weibull probability estimator
P = zeros(1, N);
for i = 1:N
    num = N - i;
    P(i) = (num + 1) / (N + 1);
end

% Weibull plot transformation
ln_sigma = log(sigma);
ln_ln_invP = log(log(1 ./ P));

% Perform linear regression to find Weibull modulus
coeffs = polyfit(ln_sigma, ln_ln_invP, 1);
m = coeffs(1); % Weibull modulus
sigma_0 = exp(-coeffs(2) / m); % Characteristic flexural strength

% Plot Weibull distribution
figure;
plot(ln_sigma, ln_ln_invP, 'bo', 'MarkerFaceColor', 'b'); % Data points
hold on;
plot(ln_sigma, polyval(coeffs, ln_sigma), 'r-', 'LineWidth', 2); % Fit line
grid on;
xlabel('ln(sigma)');
ylabel('ln(ln(1/P))');
title('Weibull Plot Problem 3');
legend('Data Points', 'Linear Fit', 'Location', 'best');

% Display results
fprintf('Weibull Modulus (m): %.2f\n', m);
fprintf('Characteristic Flexural Strength (σ₀): %.2f MPa\n', sigma_0);
