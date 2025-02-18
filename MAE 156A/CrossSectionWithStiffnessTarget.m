clear all;
close all;
clc;
format long;

%% This code provides the most optimal cross section of a material with given stiffness targets.
%Version 3

%% Changeable values

% Material properties (example values)
% E = 205e9;  % Young's Modulus for chromoly steel in Pa
% G = 80e9;   % Shear Modulus for chromoly steel in Pa

E = 70e9;   % Young's Modulus for aluminum in Pa
G = 26e9;   % Shear Modulus for aluminum in Pa

% Target flexural stiffness and torsional rigidity
EI_target = 1.23555e4;  % Target flexural stiffness in Nm^2
GJ_target = 5.636e3; % Target torsional rigidity in Nm^2

%Define constraints in inches
   lb_rec = [0.2, 0.2, 0.1];  %lower bound b, h and t
   ub_rec = [3, 3, 0.5];   %upper bound b, h and t

   lb_cir = [0.2, 0.1];  %lower bound r and t
   ub_cir = [3, 0.5];   %upper bound r and t

   lb_eli = [0.2, 0.2, 0.3];  %lower bound a (major), b (minor) and t
   ub_eli = [3, 3, 0.5];   %upper bound a (major), b (minor) and t

%% Rest of code
% Optimize rectangular cross-section
[b_rect, h_rect, t_rect, error_rect_EI, error_rect_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G, lb_rec, ub_rec);

% Optimize circular cross-section
[r_circ, t_circ, error_circ_EI, error_circ_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G, lb_cir, ub_cir);

% Optimize half-elliptical cross-section
[a_ellip, b_ellip, t_ellip, error_ellip_EI, error_ellip_GJ] = optimize_half_elliptical_cross_section(EI_target, GJ_target, E, G, lb_eli, ub_eli);

% Convert results to inches (1 meter = 39.3701 inches)
b_rect_inch = b_rect / 0.0254;
h_rect_inch = h_rect / 0.0254;
t_rect_inch = t_rect / 0.0254;
r_circ_inch = r_circ / 0.0254;
t_circ_inch = t_circ / 0.0254;
a_ellip_inch = a_ellip / 0.0254;
b_ellip_inch = b_ellip / 0.0254;
t_ellip_inch = t_ellip / 0.0254;

% Print results
fprintf('--- Hollow Rectangular Cross-Section ---\n');
fprintf('Width (b): %.4f inches\n', b_rect_inch);
fprintf('Height (h): %.4f inches\n', h_rect_inch);
fprintf('Thickness (t): %.4f inches\n', t_rect_inch);
fprintf('Error in EI: %.4f%%\n', error_rect_EI);
fprintf('Error in GJ: %.4f%%\n', error_rect_GJ);

fprintf('\n--- Hollow Circular Cross-Section ---\n');
fprintf('Radius (r): %.4f inches\n', r_circ_inch);
fprintf('Thickness (t): %.4f inches\n', t_circ_inch);
fprintf('Error in EI: %.4f%%\n', error_circ_EI);
fprintf('Error in GJ: %.4f%%\n', error_circ_GJ);

fprintf('\n--- Half-Elliptical Cross-Section ---\n');
fprintf('Major Axis (a): %.4f inches\n', a_ellip_inch);
fprintf('Minor Axis (b): %.4f inches\n', b_ellip_inch);
fprintf('Thickness (t): %.4f inches\n', t_ellip_inch);
fprintf('Error in EI: %.4f%%\n', error_ellip_EI);
fprintf('Error in GJ: %.4f%%\n', error_ellip_GJ);

%% Function to optimize the rectangular cross-section
function [b, h, t, error_EI, error_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G, lb_rec, ub_rec);
    % Error function to minimize
    errorFunc_rect = @(x) abs(log(calculate_EI_rect(x(1), x(2), x(3), E) / EI_target)) + ...
                     abs(log(calculate_GJ_rect(x(1), x(2), x(3), G) / GJ_target));

    
    % Initial guess for b, h, and t (in meters)
    initial_guess_rect = [0.05, 0.05, 0.001]; % b, h, t

    % Bounds
    lb = lb_rec * 0.0254;  % Small nonzero values
    ub = ub_rec * 0.0254;  % Larger values

    % Optimization options
    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');

    % Use fmincon for bounded optimization
    optimal_params_rect = fmincon(errorFunc_rect, initial_guess_rect, [], [], [], [], lb, ub, [], options);

    % Extract optimized values
    b = optimal_params_rect(1);
    h = optimal_params_rect(2);
    t = optimal_params_rect(3);

    % Calculate the actual EI and GJ for the optimized rectangle
    EI_actual = calculate_EI_rect(b, h, t, E);
    GJ_actual = calculate_GJ_rect(b, h, t, G);

    % Calculate the percentage errors in EI and GJ
    error_EI = ((EI_actual - EI_target) / EI_target) * 100;
    error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
end

%% Function to optimize the circular cross-section
function [r, t, error_EI, error_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G, lb_cir, ub_cir);
    errorFunc_circ = @(x) (calculate_EI_circ(x(1), x(2), E) - EI_target)^2 + ...
                           (calculate_GJ_circ(x(1), x(2), G) - GJ_target)^2;
    
    initial_guess_circ = [0.05, 0.001]; % Initial guesses: r, t

    % Bounds: ensure r > t, and practical limits
    lb = lb_cir* 0.0254;  % Smallest reasonable values
    ub = ub_cir * 0.0254;   % r ≤ 3 inches, t ≤ 0.01 meters

    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
    optimal_params_circ = fmincon(errorFunc_circ, initial_guess_circ, [], [], [], [], lb, ub, [], options);
    
    r = optimal_params_circ(1);
    t = optimal_params_circ(2);

    EI_actual = calculate_EI_circ(r, t, E);
    GJ_actual = calculate_GJ_circ(r, t, G);
    
    error_EI = ((EI_actual - EI_target) / EI_target) * 100;
    error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
end

%% Function to optimize the half-elliptical cross-section
function [a, b, t, error_EI, error_GJ] = optimize_half_elliptical_cross_section(EI_target, GJ_target, E, G, lb_eli, ub_eli);
    errorFunc_ellip = @(x) (calculate_EI_ellip(x(1), x(2), x(3), E) - EI_target)^2 + ...
                           (calculate_GJ_ellip(x(1), x(2), x(3), G) - GJ_target)^2;
    
    initial_guess_ellip = [0.05, 0.025, 0.001]; % Initial guess: a, b, t

    % Bounds: a (major axis), b (minor axis) ≤ 3 inches, thickness reasonable
    lb = lb_eli * 0.0254;  % Smallest values
    ub = ub_eli * 0.0254;  % largest values

    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
    optimal_params_ellip = fmincon(errorFunc_ellip, initial_guess_ellip, [], [], [], [], lb, ub, [], options);
    
    a = optimal_params_ellip(1);
    b = optimal_params_ellip(2);
    t = optimal_params_ellip(3);

    EI_actual = calculate_EI_ellip(a, b, t, E);
    GJ_actual = calculate_GJ_ellip(a, b, t, G);
    
    error_EI = ((EI_actual - EI_target) / EI_target) * 100;
    error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
end


%% Function to calculate EI for a rectangular section
function EI = calculate_EI_rect(b, h, t, E)
    b_inner = b - 2 * t; 
    h_inner = h - 2 * t;
    if b_inner <= 0 || h_inner <= 0
        EI = Inf; % Penalize infeasible values
    else
        EI = (E / 12) * (b * h^3 - b_inner * h_inner^3);
    end
end

%% Function to calculate GJ for a rectangular section
function GJ = calculate_GJ_rect(b, h, t, G)
    b_inner = b - 2 * t;  
    h_inner = h - 2 * t;
    if b_inner <= 0 || h_inner <= 0
        GJ = Inf; % Penalize infeasible values
    else
        GJ = G * ( (1/12*b*h*(b^2+h^2)) - (1/12*b_inner*h_inner*(b_inner^2+h_inner^2)) );
    end
end

%% Function to calculate EI for a circular section
function EI = calculate_EI_circ(r, t, E)
    r_inner = r - t;
    if r_inner <= 0
        EI = Inf; % Penalize infeasible values
    else
        EI = (E * pi / 4) * (r^4 - r_inner^4);
    end
end

%% Function to calculate GJ for a circular section
function GJ = calculate_GJ_circ(r, t, G)
    r_inner = r - t;
    if r_inner <= 0
        GJ = Inf; % Penalize infeasible values
    else
        GJ = (G * pi / 2) * (r^4 - r_inner^4);
    end
end


%% Function to calculate EI for a half-elliptical section
function EI = calculate_EI_ellip(a, b, t, E)
    a_inner = a - t;
    b_inner = b - t;
    if a_inner <= 0 || b_inner <= 0
        EI = Inf;
    else
        % Second moment of area for a half-ellipse about the flat base
        I_outer = (pi / 8 - 8/(9*pi)) * (a * b^3);
        I_inner = (pi / 8 - 8/(9*pi)) * (a_inner * b_inner^3);
        EI = E * (I_outer - I_inner);
    end
end

%% Function to calculate GJ for a half-elliptical section
function GJ = calculate_GJ_ellip(a, b, t, G)
    a_inner = a - t;
    b_inner = b - t;
    if a_inner <= 0 || b_inner <= 0
        GJ = Inf;
    else
        % Approximate polar moment of inertia for a half-ellipse
        J_outer = (pi / 8 - 8/(9*pi)) * (a * b^3 + a^3 * b);
        J_inner = (pi / 8 - 8/(9*pi)) * (a_inner * b_inner^3 + a_inner^3 * b_inner);
        GJ = G * (J_outer - J_inner);
    end
end





% %Version 1
% 
% % Material properties (example values)
% E = 210e9;  % Young's Modulus for steel in Pa
% G = 80e9;   % Shear Modulus for steel in Pa
% 
% % Target flexural stiffness and torsional rigidity (example values)
% EI_target = 7.688e3; % Target flexural stiffness in Nm^2
% GJ_target = 1.9579e4; % Target torsional rigidity in Nm^2
% 
% % Defined Thickness
% t = 0.1 * 0.0254; % 0.1 inches to meters
% 
% % Find the optimized cross-sections for rectangle and circle
% [b_rect, h_rect, error_rect_EI, error_rect_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G, t);
% [r_circ, error_circ_EI, error_circ_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G, t);
% 
% % Convert the results from meters to inches (1 meter = 39.3701 inches)
% b_rect_inch = b_rect / 0.0254;
% h_rect_inch = h_rect / 0.0254;
% r_circ_inch = r_circ / 0.0254;
% 
% % Print the results in inches
% fprintf('--- Hollow Rectangular Cross-Section ---\n');
% fprintf('Optimized Width (b): %.4f inches\n', b_rect_inch);
% fprintf('Optimized Height (h): %.4f inches\n', h_rect_inch);
% fprintf('Error in Flexural Stiffness (EI): %.4f%%\n', error_rect_EI);
% fprintf('Error in Torsional Rigidity (GJ): %.4f%%\n', error_rect_GJ);
% 
% fprintf('\n--- Hollow Circular Cross-Section ---\n');
% fprintf('Optimized Radius (r): %.4f inches\n', r_circ_inch);
% fprintf('Error in Flexural Stiffness (EI): %.4f%%\n', error_circ_EI);
% fprintf('Error in Torsional Rigidity (GJ): %.4f%%\n', error_circ_GJ);
% 
% %% Function to optimize the rectangular cross-section
% function [b, h, error_EI, error_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G, t)
%     errorFunc_rect = @(x) (calculate_EI_rect(x(1), x(2), E, t) - EI_target)^2 + ...
%                            (calculate_GJ_rect(x(1), x(2), G, t) - GJ_target)^2;
%     
%     initial_guess_rect = [0.1, 0.1]; % Reasonable initial guess in meters
%     
%     options = optimset('TolX', 1e-6, 'TolFun', 1e-6);
%     optimal_params_rect = fminsearch(errorFunc_rect, initial_guess_rect, options);
%     
%     b = optimal_params_rect(1);
%     h = optimal_params_rect(2);
%     
%     EI_actual = calculate_EI_rect(b, h, E, t);
%     GJ_actual = calculate_GJ_rect(b, h, G, t);
%     
%     error_EI = ((EI_actual - EI_target) / EI_target) * 100;
%     error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
% end
% 
% %% Function to optimize the circular cross-section
% function [r, error_EI, error_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G, t)
%     errorFunc_circ = @(r) (calculate_EI_circ(r, E, t) - EI_target)^2 + ...
%                            (calculate_GJ_circ(r, G, t) - GJ_target)^2;
%     
%     initial_guess_circ = 0.1; % Reasonable initial guess in meters
%     
%     options = optimset('TolX', 1e-6, 'TolFun', 1e-6);
%     r_optimal = fminsearch(errorFunc_circ, initial_guess_circ, options);
%     
%     r = r_optimal;
%     
%     EI_actual = calculate_EI_circ(r, E, t);
%     GJ_actual = calculate_GJ_circ(r, G, t);
%     
%     error_EI = ((EI_actual - EI_target) / EI_target) * 100;
%     error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
% end
% 
% %% Function to calculate EI for a rectangular section
% function EI = calculate_EI_rect(b, h, E, t)
%     b_inner = b - 2 * t; 
%     h_inner = h - 2 * t;
%     if b_inner <= 0 || h_inner <= 0
%         EI = Inf; % Penalize infeasible values
%     else
%         EI = (E / 12) * (b * h^3 - b_inner * h_inner^3);
%     end
% end
% 
% %% Function to calculate GJ for a rectangular section
% function GJ = calculate_GJ_rect(b, h, G, t)
%     b_inner = b - 2 * t;  
%     h_inner = h - 2 * t;
%     if b_inner <= 0 || h_inner <= 0
%         GJ = Inf; % Penalize infeasible values
%     else
%         % Approximate torsional constant
%         GJ = (4 * b * h^3 * G / 3) * ((1 - (b_inner * h_inner^3) / (b * h^3)) / (1 + h / b));
%     end
% end
% 
% %% Function to calculate EI for a circular section
% function EI = calculate_EI_circ(r, E, t)
%     r_inner = r - t;
%     if r_inner <= 0
%         EI = Inf; % Penalize infeasible values
%     else
%         EI = (E * pi / 4) * (r^4 - r_inner^4);
%     end
% end
% 
% %% Function to calculate GJ for a circular section
% function GJ = calculate_GJ_circ(r, G, t)
%     r_inner = r - t;
%     if r_inner <= 0
%         GJ = Inf; % Penalize infeasible values
%     else
%         GJ = (G * pi / 2) * (r^4 - r_inner^4);
%     end
% end

% %Version 2
% % Material properties (example values)
% E = 210e9;  % Young's Modulus for steel in Pa
% G = 80e9;   % Shear Modulus for steel in Pa
% 
% % Target flexural stiffness and torsional rigidity
% EI_target = 7.688e3;  % Target flexural stiffness in Nm^2
% GJ_target = 1.9579e4; % Target torsional rigidity in Nm^2
% 
% % Find the optimized cross-sections for rectangle and circle
% [b_rect, h_rect, t_rect, error_rect_EI, error_rect_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G);
% [r_circ, t_circ, error_circ_EI, error_circ_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G);
% 
% % Convert the results from meters to inches (1 meter = 39.3701 inches)
% b_rect_inch = b_rect / 0.0254;
% h_rect_inch = h_rect / 0.0254;
% t_rect_inch = t_rect / 0.0254;
% r_circ_inch = r_circ / 0.0254;
% t_circ_inch = t_circ / 0.0254;
% 
% % Print the results in inches
% fprintf('--- Hollow Rectangular Cross-Section ---\n');
% fprintf('Optimized Width (b): %.4f inches\n', b_rect_inch);
% fprintf('Optimized Height (h): %.4f inches\n', h_rect_inch);
% fprintf('Optimized Thickness (t): %.4f inches\n', t_rect_inch);
% fprintf('Error in Flexural Stiffness (EI): %.4f%%\n', error_rect_EI);
% fprintf('Error in Torsional Rigidity (GJ): %.4f%%\n', error_rect_GJ);
% 
% fprintf('\n--- Hollow Circular Cross-Section ---\n');
% fprintf('Optimized Radius (r): %.4f inches\n', r_circ_inch);
% fprintf('Optimized Thickness (t): %.4f inches\n', t_circ_inch);
% fprintf('Error in Flexural Stiffness (EI): %.4f%%\n', error_circ_EI);
% fprintf('Error in Torsional Rigidity (GJ): %.4f%%\n', error_circ_GJ);
% 
% %% Function to optimize the rectangular cross-section
% function [b, h, t, error_EI, error_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G)
%     % Error function to minimize
%     errorFunc_rect = @(x) (calculate_EI_rect(x(1), x(2), x(3), E) - EI_target)^2 + ...
%                            (calculate_GJ_rect(x(1), x(2), x(3), G) - GJ_target)^2;
%     
%     % Initial guess for b, h, and t (in meters)
%     initial_guess_rect = [0.05, 0.05, 0.001]; % b, h, t
% 
%     % Bounds: ensure 0 < t < min(b,h)/2, and b, h ≤ 3 inches (0.0762 meters)
%     lb = [0.001, 0.001, 0.00254];  % Small nonzero values
%     ub = [0.0762, 0.0762, 0.01];  % b, h ≤ 3 inches, t ≤ 0.01 meters (0.39 inches)
% 
%     % Optimization options
%     options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
% 
%     % Use fmincon for bounded optimization
%     optimal_params_rect = fmincon(errorFunc_rect, initial_guess_rect, [], [], [], [], lb, ub, [], options);
% 
%     % Extract optimized values
%     b = optimal_params_rect(1);
%     h = optimal_params_rect(2);
%     t = optimal_params_rect(3);
% 
%     % Calculate the actual EI and GJ for the optimized rectangle
%     EI_actual = calculate_EI_rect(b, h, t, E);
%     GJ_actual = calculate_GJ_rect(b, h, t, G);
% 
%     % Calculate the percentage errors in EI and GJ
%     error_EI = ((EI_actual - EI_target) / EI_target) * 100;
%     error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
% end
% 
% %% Function to optimize the circular cross-section
% function [r, t, error_EI, error_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G)
%     errorFunc_circ = @(x) (calculate_EI_circ(x(1), x(2), E) - EI_target)^2 + ...
%                            (calculate_GJ_circ(x(1), x(2), G) - GJ_target)^2;
%     
%     initial_guess_circ = [0.05, 0.001]; % Initial guesses: r, t
% 
%     % Bounds: ensure r > t, and practical limits
%     lb = [0.005, 0.00254];  % Smallest reasonable values
%     ub = [0.0762, 0.01];   % r ≤ 3 inches, t ≤ 0.01 meters
% 
%     options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
%     optimal_params_circ = fmincon(errorFunc_circ, initial_guess_circ, [], [], [], [], lb, ub, [], options);
%     
%     r = optimal_params_circ(1);
%     t = optimal_params_circ(2);
% 
%     EI_actual = calculate_EI_circ(r, t, E);
%     GJ_actual = calculate_GJ_circ(r, t, G);
%     
%     error_EI = ((EI_actual - EI_target) / EI_target) * 100;
%     error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
% end
% 
% %% Function to calculate EI for a rectangular section
% function EI = calculate_EI_rect(b, h, t, E)
%     b_inner = b - 2 * t; 
%     h_inner = h - 2 * t;
%     if b_inner <= 0 || h_inner <= 0
%         EI = Inf; % Penalize infeasible values
%     else
%         EI = (E / 12) * (b * h^3 - b_inner * h_inner^3);
%     end
% end
% 
% %% Function to calculate GJ for a rectangular section
% function GJ = calculate_GJ_rect(b, h, t, G)
%     b_inner = b - 2 * t;  
%     h_inner = h - 2 * t;
%     if b_inner <= 0 || h_inner <= 0
%         GJ = Inf; % Penalize infeasible values
%     else
%         GJ = (4 * b * h^3 * G / 3) * ((1 - (b_inner * h_inner^3) / (b * h^3)) / (1 + h / b));
%     end
% end
% 
% %% Function to calculate EI for a circular section
% function EI = calculate_EI_circ(r, t, E)
%     r_inner = r - t;
%     if r_inner <= 0
%         EI = Inf; % Penalize infeasible values
%     else
%         EI = (E * pi / 4) * (r^4 - r_inner^4);
%     end
% end
% 
% %% Function to calculate GJ for a circular section
% function GJ = calculate_GJ_circ(r, t, G)
%     r_inner = r - t;
%     if r_inner <= 0
%         GJ = Inf; % Penalize infeasible values
%     else
%         GJ = (G * pi / 2) * (r^4 - r_inner^4);
%     end
% end
% 

