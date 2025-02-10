%Temporary Scripts
clc; clear all; close all;

%% Problem 3
% Given data
time = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]; % Time in minutes
strain = [0.000, 0.025, 0.043, 0.065, 0.078, 0.092, 0.109, 0.120, 0.135, 0.153, 0.172, 0.193, 0.218, 0.255, 0.307, 0.368]; % Strain

% Plot strain vs. time
figure;
plot(time, strain, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Strain', 'FontSize', 12);
title('Creep Strain vs. Time for Aluminum Alloy at 400°C', 'FontSize', 14);
grid on;

% Determine the steady-state creep rate (secondary region)
% Select points in the steady-state region (approximately linear portion)
steady_time = time(8:12); % Adjust indices if needed
steady_strain = strain(8:12);

% Linear regression
coeffs = polyfit(steady_time, steady_strain, 1);
steady_state_creep_rate = coeffs(1); % Slope of strain vs. time

% Display results
disp(['Steady-state creep rate: ', num2str(steady_state_creep_rate), ' per minute']);

% Add trendline to plot
hold on;
plot(steady_time, polyval(coeffs, steady_time), 'r--', 'LineWidth', 2);
legend('Creep Data', 'Steady-State Fit', 'Location', 'NorthWest');


%% Problem 5
% Given data
T_C = { [400, 600, 800, 1000], ...  % High-lead glass
        [600, 800, 1000, 1200, 1400, 1600], ... % Soda-lime glass
        [1200, 1400, 1600] }; % Fused silica

viscosity = { [1e14, 3.16227766e7, 3.16227766e4, 1000], ... % High-lead glass
              [1e11, 3.16227766e6, 3.16227766e4, 3162.27766, 501.1872336, 125.8925412], ... % Soda-lime glass
              [3.16228e12, 1e10, 1.258925412e8] }; % Fused silica
          
R = 8.314; % Gas constant in J/(mol*K)

figure;
colors = ['r', 'g', 'b'];
labels = {'High-lead glass', 'Soda-lime glass', 'Fused silica'};
Ea_values = zeros(1,3); % To store activation energies

for i = 1:3
    T_K = T_C{i} + 273.15; % Convert to Kelvin
    inv_T = 1 ./ T_K; % Compute 1/T
    ln_visc = log(viscosity{i}); % Natural log of viscosity
    
    % Perform linear regression
    coeffs = polyfit(inv_T, ln_visc, 1);
    slope = coeffs(1);
    
    % Compute activation energy Q
    Q = slope * R;
    Ea_values(i) = Q;
    
    % Plot results
    subplot(1,3,i);
    plot(inv_T, ln_visc, 'o', 'Color', colors(i), 'MarkerFaceColor', colors(i));
    hold on;
    plot(inv_T, polyval(coeffs, inv_T), 'Color', colors(i));
    xlabel('1/T (1/K)');
    ylabel('ln(viscosity)');
    title(labels{i});
    grid on;
end

disp('Activation energies (J/mol):');
disp(table(labels', Ea_values', 'VariableNames', {'Glass Type', 'Activation Energy (J/mol)'}));





















% 
% A = [1/2 1/2 -1/sqrt(2); -1/sqrt(2) 1/sqrt(2) 0; 1/2 1/2 1/sqrt(2)]
% 
% B = [1-sqrt(2) 1 0; 1 1+sqrt(2) 0; 0 0 -2]
% 
% C = [1/2 -1/sqrt(2) 1/2; 1/2 1/sqrt(2) 1/2; -1/sqrt(2) 0 1/sqrt(2)]
% 
% 
% F = tf([1], [1 10 24 0])
% % allmargin(tf)
% rlocus(F)
% 
% z=0.504;
% sgrid(z,0);
% 
% 
% % Data
% x_values = [45, 41, 37, 33, 29, 25, 21, 17, 13, 9];  % distance in cm
% y_values = [2500, 2270, 2030, 1800, 1590, 1380, 1110, 890, 730, 550];  % duration in microseconds
% 
% % Plot the data points
% figure;
% scatter(x_values, y_values, 'b', 'filled');
% hold on;
% 
% % Fit a line to the data
% p = polyfit(x_values, y_values, 1);
% y_fit = polyval(p, x_values);
% 
% % Plot the line of best fit
% plot(x_values, y_fit, 'r-', 'LineWidth', 1.5);
% 
% % Labels and title
% xlabel('Distance (cm)');
% ylabel('Duration (μs)');
% 
% grid on;
% 
% 
% 
% % Speed of sound calculation
% speed_of_sound_exp = 1 / p(1) * 10000;  % converting from cm/μs to m/s
% speed_of_sound_theoretical = 343;  % theoretical speed of sound in air at 20°C in m/s
% 
% % Comparison
% fprintf('Experimental Speed of Sound: %.2f m/s\n', speed_of_sound_exp);
% fprintf('Theoretical Speed of Sound: %.2f m/s\n', speed_of_sound_theoretical);
% fprintf('Difference: %.2f m/s\n', abs(speed_of_sound_exp - speed_of_sound_theoretical));
% 
% hold off;

% % Define symbolic variables
% syms P b L E w h real
% 
% % Moment of inertia
% I = (1/12) * w * h^3;
% 
% % Deflection formula
% Vmax = (P * b * (L^2 - b^2)^(3/2)) / (9 * sqrt(3) * E * I * L);
% 
% % Derivatives
% dVmax_dP = diff(Vmax, P);
% dVmax_db = diff(Vmax, b);
% dVmax_dL = diff(Vmax, L);
% dVmax_dE = diff(Vmax, E);
% dVmax_dw = diff(Vmax, w);
% dVmax_dh = diff(Vmax, h);
% 
% % Display the symbolic derivatives
% disp('dVmax/dP:');
% pretty(dVmax_dP)
% 
% disp('dVmax/db:');
% pretty(dVmax_db)
% 
% disp('dVmax/dL:');
% pretty(dVmax_dL)
% 
% disp('dVmax/dE:');
% pretty(dVmax_dE)
% 
% disp('dVmax/dw:');
% pretty(dVmax_dw)
% 
% disp('dVmax/dh:');
% pretty(dVmax_dh)
% 
% % Given values with uncertainties
% L_val = 1;         % m
% delta_L = 0.01;    % m
% w_val = 0.05;      % m (5 cm)
% delta_w = 0.002;   % m (0.2 cm)
% h_val = 0.03;      % m (3 cm)
% delta_h = 0.002;   % m (0.2 cm)
% E_val = 107e9;     % Pa (107 GPa)
% delta_E = 5e9;     % Pa (5 GPa)
% P_val = 6000;      % N (6 kN)
% delta_P = 2000;    % N (2 kN)
% b_val = 0.3;       % m
% delta_b = 0.1;     % m
% 
% % Recalculate I with nominal values
% I_val = (1/12) * w_val * h_val^3;
% 
% % Nominal value of maximum deflection Vmax
% Vmax_nominal = (P_val * b_val * (L_val^2 - b_val^2)^(3/2)) / (9 * sqrt(3) * E_val * I_val * L_val);
% 
% % Evaluate each derivative numerically
% dVmax_dP_num = double(subs(dVmax_dP, {P, b, L, E, w, h}, {P_val, b_val, L_val, E_val, w_val, h_val}));
% dVmax_db_num = double(subs(dVmax_db, {P, b, L, E, w, h}, {P_val, b_val, L_val, E_val, w_val, h_val}));
% dVmax_dL_num = double(subs(dVmax_dL, {P, b, L, E, w, h}, {P_val, b_val, L_val, E_val, w_val, h_val}));
% dVmax_dE_num = double(subs(dVmax_dE, {P, b, L, E, w, h}, {P_val, b_val, L_val, E_val, w_val, h_val}));
% dVmax_dw_num = double(subs(dVmax_dw, {P, b, L, E, w, h}, {P_val, b_val, L_val, E_val, w_val, h_val}));
% dVmax_dh_num = double(subs(dVmax_dh, {P, b, L, E, w, h}, {P_val, b_val, L_val, E_val, w_val, h_val}));
% 
% % Print numerical results
% fprintf('dVmax/dP = %.6e\n', dVmax_dP_num);
% fprintf('dVmax/db = %.6e\n', dVmax_db_num);
% fprintf('dVmax/dL = %.6e\n', dVmax_dL_num);
% fprintf('dVmax/dE = %.6e\n', dVmax_dE_num);
% fprintf('dVmax/dw = %.6e\n', dVmax_dw_num);
% fprintf('dVmax/dh = %.6e\n', dVmax_dh_num);
% 
% % Error propagation calculation
% delta_Vmax = sqrt((dVmax_dP_num * delta_P)^2 + ...
%                   (dVmax_db_num * delta_b)^2 + ...
%                   (dVmax_dL_num * delta_L)^2 + ...
%                   (dVmax_dE_num * delta_E)^2 + ...
%                   (dVmax_dw_num * delta_w)^2 + ...
%                   (dVmax_dh_num * delta_h)^2);
% 
% % Print final result
% fprintf('Nominal Maximum Deflection: %.6e m\n', Vmax_nominal);
% fprintf('Uncertainty in Maximum Deflection: %.6e m\n', delta_Vmax);
% 
% 
% 
% MATLAB Script to Symbolically and Numerically Solve for Maximum Concentration Time
% 
% % % Clear previous variables and close figures
% clear; clc; close all;
% 
% % Define symbolic variables for distance (x), time (t), and constants
% syms x t M D U L A
% 
% % Define the Gaussian concentration function C(x, t)
% C = (M / (A*sqrt(4 * pi * D * t))) * exp(-(x - U * t)^2 / (4 * D * t));
% 
% % Compute the derivative of C with respect to t
% dC_dt = diff(C, t);
% dC_dt2 = diff(diff(C, x), x).*D;
% % Display the derivative
% disp('Time derivative of concentration, dC/dt:');
% pretty(dC_dt);
% 
% % Solve dC/dt = 0 to find the time t that maximizes C at x = L
% t_max_sol = solve(subs(dC_dt, x, L) == 0, t);
% t_max_sol2 = solve(subs(dC_dt2, x, L) == 0, t);
% 
% % Display the symbolic solution for the maximum concentration time
% disp('Symbolic solution for the time at which concentration is maximized:');
% disp(t_max_sol);
% 
% % Substitute numerical values for the parameters
% M_val = 10;        % kg of contaminant
% D_val = 0.1;       % m^2/s, turbulent diffusivity
% U_val = 0.007;     % m/s, flow velocity
% L_val = -100;       % m, distance to hatchery
% A_val = 20;        %m^2, cross sectional area
% 
% % Substitute values into the symbolic solution and evaluate numerically
% t_max_num = double(subs(t_max_sol, [M, D, U, L, A], [M_val, D_val, U_val, L_val, A_val]));
% 
% % Display the numerical result for the time of maximum concentration
% disp('Numerical value for the time at which concentration is maximized (s):');
% disp(t_max_num(1));
% 
% % Compute the maximum concentration at the time t = t_max_num
% C_max_num = double(subs(C, {x, t, M, D, U, L, A}, {L_val, t_max_num(1), M_val, D_val, U_val, L_val, A_val}));
% 
% % Display the numerical result for the maximum concentration
% disp('Numerical value for the maximum concentration (kg/m^3):');
% disp(C_max_num);



% % Constants
% R = 8.314; % Universal gas constant in J/(K*mol)
% Ea_nom = 40000; % Nominal activation energy in J/mol
% Ea_unc = 10000; % Uncertainty in activation energy in J/mol
% T_nom = 289; % Nominal temperature in K
% T_unc = 5; % Uncertainty in temperature in K
% beta_nom = 0.8; % Nominal fudge factor
% beta_unc = 0.1; % Uncertainty in fudge factor
% k0_nom = 1.5e6; % Nominal attempt frequency in s^-1
% k0_unc = 0.5e6; % Uncertainty in attempt frequency in s^-1
% 
% % Calculate nominal k
% k_nom = k0_nom * exp(-(Ea_nom / (R * T_nom))^beta_nom);
% 
% % Partial derivatives for uncertainty propagation
% % Derivative with respect to Ea
% dEak = k0_nom * exp(-(Ea_nom / (R * T_nom))^beta_nom) * (-beta_nom * (Ea_nom / (R * T_nom))^(beta_nom - 1)) / (R * T_nom);
% % Derivative with respect to T
% dTk = k0_nom * exp(-(Ea_nom / (R * T_nom))^beta_nom) * (beta_nom * Ea_nom / (R * T_nom^2)) * (Ea_nom / (R * T_nom))^(beta_nom - 1);
% % Derivative with respect to beta
% dbetak = -k0_nom * exp(-(Ea_nom / (R * T_nom))^beta_nom) * log(Ea_nom / (R * T_nom)) * (Ea_nom / (R * T_nom))^beta_nom;
% % Derivative with respect to k0
% dk0k = exp(-(Ea_nom / (R * T_nom))^beta_nom);
% 
% % Calculate uncertainty in k
% delta_k = sqrt((dEak * Ea_unc)^2 + (dTk * T_unc)^2 + (dbetak * beta_unc)^2 + (dk0k * k0_unc)^2);
% 
% % Display results
% fprintf('Nominal value of k: %.3e s^-1\n', k_nom);
% fprintf('Uncertainty in k: ± %.3e s^-1\n', delta_k);

% % Parameters
% L = 5;                % Total rise height in mm
% rise_time = 0.02;     % Duration of rise phase in seconds
% dwell_time = 0.1;     % Duration of dwell phase in seconds
% fall_time = 0.02;     % Duration of fall phase in seconds
% total_time = rise_time + dwell_time + fall_time; % Total time for one cycle
% 
% % Constant angular velocity
% omega = 2 * pi / total_time; % Angular speed for one complete cycle in rad/s
% 
% % Define time intervals for each phase
% t_rise = linspace(0, rise_time, 100);                  % Time for rise phase
% t_dwell = linspace(rise_time, rise_time + dwell_time, 100); % Time for dwell phase
% t_fall = linspace(rise_time + dwell_time, total_time, 100); % Time for fall phase
% 
% % Calculate angular position (theta) for each phase
% theta_rise = omega * t_rise;                          % Angular position during rise
% theta_dwell = omega * t_dwell;                        % Angular position during dwell
% 
% % Set theta_fall equal to theta_rise to ensure symmetry
% theta_fall = theta_rise; 
% 
% % Calculate displacement, velocity, and acceleration for each phase
% % Using the provided cycloidal rise profile y = L[theta/beta - (1/(2*pi)) * sin(2*pi*theta/beta)]
% 
% beta = omega * rise_time; % Total angular displacement for the rise phase in radians
% 
% % Rise phase
% y_rise = L * (theta_rise / beta - (1 / (2 * pi)) * sin(2 * pi * theta_rise / beta));
% v_rise = L / beta * (1 - cos(2 * pi * theta_rise / beta)) * omega; % dy/dtheta * omega
% a_rise = L * (2 * pi / beta^2) * sin(2 * pi * theta_rise / beta) * omega^2; % d2y/dtheta2 * omega^2
% 
% % Dwell phase (constant displacement, zero velocity and acceleration)
% y_dwell = L * ones(size(t_dwell)); % Displacement remains constant during dwell
% v_dwell = zeros(size(t_dwell));    % Velocity is zero during dwell
% a_dwell = zeros(size(t_dwell));    % Acceleration is zero during dwell
% 
% % Fall phase (mirror the rise phase to return to zero displacement)
% y_fall = L * (1 - (theta_fall / beta - (1 / (2 * pi)) * sin(2 * pi * theta_fall / beta)));
% v_fall = -L / beta * (1 - cos(2 * pi * theta_fall / beta)) * omega; % Negative of rise velocity
% a_fall = -L * (2 * pi / beta^2) * sin(2 * pi * theta_fall / beta) * omega^2; % Negative of rise acceleration
% 
% % Concatenate results for full cycle
% time = [t_rise, t_dwell, t_fall];
% displacement = [y_rise, y_dwell, y_fall];
% velocity = [v_rise, v_dwell, v_fall];
% acceleration = [a_rise, a_dwell, a_fall];
% 
% % Plot the results
% figure;
% subplot(3,1,1);
% plot(time, displacement, 'b');
% title('Displacement Profile');
% ylabel('Displacement y (mm)');
% 
% subplot(3,1,2);
% plot(time, velocity, 'g');
% title('Velocity Profile');
% ylabel('Velocity v (mm/s)');
% 
% subplot(3,1,3);
% plot(time, acceleration, 'r');
% title('Acceleration Profile');
% ylabel('Acceleration a (mm/s^2)');
% xlabel('Time (s)');
