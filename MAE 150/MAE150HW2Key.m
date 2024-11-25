clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1%
% (a) Load the tabulated data from the file TP.dat and create a scatter plot.
data = readmatrix('TP.dat');  % Reads the data from TP.dat, ignoring the header
T = data(:, 1);  % Temperature (x-axis)
P = data(:, 2);  % Power (y-axis)

figure;
scatter(T, P, 'r', 'filled'); % Scatter plot of raw data
xlabel('Temperature (K)');
ylabel('Power (W)');
grid on;
hold on; % Keep plot open for overlaying fits

% (b) Determine the regression coefficients for different fits (Linear, Quadratic, Cubic, Quartic)
n = length(T);

% For each fit type, construct matrix A for normal equation A*a = b
% Linear fit: a1*T + a0
A_linear = [T ones(n,1)];
a_linear = A_linear \ P;  % Solving for coefficients

% Quadratic fit: a2*T^2 + a1*T + a0
A_quadratic = [T.^2 T ones(n,1)];
a_quadratic = A_quadratic \ P;

% Cubic fit: a3*T^3 + a2*T^2 + a1*T + a0
A_cubic = [T.^3 T.^2 T ones(n,1)];
a_cubic = A_cubic \ P;

% Quartic fit: a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0
A_quartic = [T.^4 T.^3 T.^2 T ones(n,1)];
a_quartic = A_quartic \ P;

% Display the regression coefficients
disp('Linear coefficients:'), disp(a_linear');
disp('Quadratic coefficients:'), disp(a_quadratic');
disp('Cubic coefficients:'), disp(a_cubic');
disp('Quartic coefficients:'), disp(a_quartic');

% (c) Plot best-fit approximations on the same figure with the raw data
T_fit = linspace(min(T), max(T), 100);  % Finer temperature range for smooth plotting

% Linear fit
P_linear_fit = a_linear(1)*T_fit + a_linear(2);
plot(T_fit, P_linear_fit, 'b-', 'LineWidth', 1.5);

% Quadratic fit
P_quadratic_fit = a_quadratic(1)*T_fit.^2 + a_quadratic(2)*T_fit + a_quadratic(3);
plot(T_fit, P_quadratic_fit, 'k-', 'LineWidth', 1.5);

% Cubic fit
P_cubic_fit = a_cubic(1)*T_fit.^3 + a_cubic(2)*T_fit.^2 + a_cubic(3)*T_fit + a_cubic(4);
plot(T_fit, P_cubic_fit, 'g-', 'LineWidth', 1.5);

% Quartic fit
P_quartic_fit = a_quartic(1)*T_fit.^4 + a_quartic(2)*T_fit.^3 + a_quartic(3)*T_fit.^2 + a_quartic(4)*T_fit + a_quartic(5);
plot(T_fit, P_quartic_fit, 'm-', 'LineWidth', 1.5);

% Add title, legend, and grid
title('My fits');
legend('Raw Data', 'Linear fit', 'Quadratic fit', 'Cubic fit', 'Quartic fit');

% (d) Calculate the root-mean-square error (RMSE) for each fit
rmse_linear = sqrt(mean((P - (a_linear(1)*T + a_linear(2))).^2));
rmse_quadratic = sqrt(mean((P - (a_quadratic(1)*T.^2 + a_quadratic(2)*T + a_quadratic(3))).^2));
rmse_cubic = sqrt(mean((P - (a_cubic(1)*T.^3 + a_cubic(2)*T.^2 + a_cubic(3)*T + a_cubic(4))).^2));
rmse_quartic = sqrt(mean((P - (a_quartic(1)*T.^4 + a_quartic(2)*T.^3 + a_quartic(3)*T.^2 + a_quartic(4)*T + a_quartic(5))).^2));

% Print RMSE values and determine the best fit
disp(['RMSE for Linear fit: ', num2str(rmse_linear)]);
disp(['RMSE for Quadratic fit: ', num2str(rmse_quadratic)]);
disp(['RMSE for Cubic fit: ', num2str(rmse_cubic)]);
disp(['RMSE for Quartic fit: ', num2str(rmse_quartic)]);

[~, best_fit] = min([rmse_linear, rmse_quadratic, rmse_cubic, rmse_quartic]);
best_fit_names = {'Linear', 'Quadratic', 'Cubic', 'Quartic'};
disp(['The best fit is: ', best_fit_names{best_fit}]);

% (e) Use MATLAB’s polyfit function for the same regression models and plot in a new figure
p_linear = polyfit(T, P, 1);
p_quadratic = polyfit(T, P, 2);
p_cubic = polyfit(T, P, 3);
p_quartic = polyfit(T, P, 4);

figure;
scatter(T, P, 'r', 'filled');  % Replot the raw data
hold on;

% Plot polyfit curves
plot(T_fit, polyval(p_linear, T_fit), 'b-', 'LineWidth', 1.5);  % Linear fit
plot(T_fit, polyval(p_quadratic, T_fit), 'k-', 'LineWidth', 1.5);  % Quadratic fit
plot(T_fit, polyval(p_cubic, T_fit), 'g-', 'LineWidth', 1.5);  % Cubic fit
plot(T_fit, polyval(p_quartic, T_fit), 'm-', 'LineWidth', 1.5);  % Quartic fit

% Add title, legend, and grid
title("MATLAB's polyfit");
legend('Raw Data', 'Linear fit', 'Quadratic fit', 'Cubic fit', 'Quartic fit');
xlabel('Temperature (K)');
ylabel('Power (W)');
grid on;


%%Problem 3

% (a) Load the data
data = load('xy.dat'); % Load data from the file
x = data(:, 1); % First column is x
y = data(:, 2); % Second column is y

% Generate scatter plot
figure;
scatter(x, y, 'b', 'filled'); % Blue dots for the scatter plot
grid on;
xlabel('Normalized Quantity x');
ylabel('Normalized Quantity y');
title('Scatter Plot of Data');
hold on;

% (b) Linearizing the model: y(x) = a * 10^(b*x) + 273
% Rearranging gives: (y - 273) = a * 10^(b*x)
% Taking log: log10(y - 273) = log10(a) + b*x

% (c) Linear regression coefficients calculation
y_linearized = log10(y - 273); % Linearized y values
% Perform linear regression
p = polyfit(x, y_linearized, 1); % p(1) = b, p(2) = log10(a)

% Calculate a and b from the coefficients
b = p(1);
log10_a = p(2);
a = 10^log10_a;

% Display computed values
fprintf('Computed values:\n');
fprintf('a = %.4f\n', a);
fprintf('b = %.4f\n', b);

% (d) Plot the best-fit approximation
x_fit = linspace(min(x), max(x), 100); % Generate points for the fit line
y_fit = a * 10.^(b * x_fit) + 273; % Calculate y values for the fit

plot(x_fit, y_fit, 'k-', 'LineWidth', 2); % Black line for best fit
legend('Data Points', 'Best Fit Line');
hold off;

