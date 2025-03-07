clear all;
close all;
clc;
format long;

%% This code provides the most optimal cross section of a material with given stiffness and rigidity targets.
%Version 3

%% Changeable values

% Material properties
E = 205e9;  % Young's Modulus for chromoly steel in Pa
G = 80e9;   % Shear Modulus for chromoly steel in Pa
% 
% E = 70e9;   % Young's Modulus for aluminum in Pa
% G = 26e9;   % Shear Modulus for aluminum in Pa

%Referenced inertia values
% I = 2*9.3e-9;  %chassis 
% J = 2*1.05e-8;  %chassis 
% I = 2*8.998e-9;  %chassis 
% J = 2*9.6477e-9;  %chassis 
% I = 2*3.53e-8;  %swignarm
% J = 2*3.91e-8;  %swignarm
I = 2*5.500e-7;  %chassis plate
J = 2*7.75e-7;  %chassis plate

% Target flexural stiffness and torsional rigidity
EI_target = 70e9*I;  % Target flexural stiffness in Nm^2
GJ_target = 26e9*J; % Target torsional rigidity in Nm^2

%Define constraints in inches

   lb_rec = [1, 2, 0.065];  %lower bound b, h and t
   ub_rec = [1, 2, 0.065];   %upper bound b, h and t

   lb_cir = [3/8, 0.12];  %lower bound r and t
   ub_cir = [3/8, 0.12];   %upper bound r and t

   lb_eli = [0.2, 0.2, 0.035];  %lower bound a (major), b (minor) and t
   ub_eli = [3, 3, 0.1];   %upper bound a (major), b (minor) and t

%Number of cross sectional supports. NOTE: This assumes the supports share
%the load equally.
supports = 2;


%% Rest of code
% Optimize rectangular cross-section
[b_rect, h_rect, t_rect, error_rect_EI, error_rect_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G, lb_rec, ub_rec, supports);

% Optimize circular cross-section
[r_circ, t_circ, error_circ_EI, error_circ_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G, lb_cir, ub_cir, supports);

% Optimize half-elliptical cross-section
[a_ellip, b_ellip, t_ellip, error_ellip_EI, error_ellip_GJ] = optimize_half_elliptical_cross_section(EI_target, GJ_target, E, G, lb_eli, ub_eli, supports);

% Convert results to inches (1 meter = 39.3701 inches)
b_rect_inch = b_rect / 0.0254;
h_rect_inch = h_rect / 0.0254;
t_rect_inch = t_rect / 0.0254;
r_circ_inch = r_circ / 0.0254;
t_circ_inch = t_circ / 0.0254;
a_ellip_inch = a_ellip / 0.0254;
b_ellip_inch = b_ellip / 0.0254;
t_ellip_inch = t_ellip / 0.0254;

% Cross-sectional area calculations (in square inches)
A_rect_inch = (b_rect_inch * h_rect_inch) - ((b_rect_inch - 2*t_rect_inch) * (h_rect_inch - 2*t_rect_inch));
A_circ_inch = pi * (r_circ_inch^2 - (r_circ_inch - t_circ_inch)^2);
A_ellip_inch = (pi / 2) * (a_ellip_inch * b_ellip_inch - (a_ellip_inch - t_ellip_inch) * (b_ellip_inch - t_ellip_inch));

% Print results
fprintf('--- Hollow Rectangular Cross-Section ---\n');
fprintf('Width (b): %.4f inches\n', b_rect_inch);
fprintf('Height (h): %.4f inches\n', h_rect_inch);
fprintf('Thickness (t): %.4f inches\n', t_rect_inch);
fprintf('Cross-sectional Area: %.4f square inches\n', A_rect_inch);
fprintf('Error in EI: %.4f%%\n', error_rect_EI);
fprintf('Error in GJ: %.4f%%\n', error_rect_GJ);

fprintf('\n--- Hollow Circular Cross-Section ---\n');
fprintf('Radius (r): %.4f inches\n', r_circ_inch);
fprintf('Thickness (t): %.4f inches\n', t_circ_inch);
fprintf('Cross-sectional Area: %.4f square inches\n', A_circ_inch);
fprintf('Error in EI: %.4f%%\n', error_circ_EI);
fprintf('Error in GJ: %.4f%%\n', error_circ_GJ);

fprintf('\n--- Half-Elliptical Cross-Section ---\n');
fprintf('Major Axis (a): %.4f inches\n', a_ellip_inch);
fprintf('Minor Axis (b): %.4f inches\n', b_ellip_inch);
fprintf('Thickness (t): %.4f inches\n', t_ellip_inch);
fprintf('Cross-sectional Area: %.4f square inches\n', A_ellip_inch);
fprintf('Error in EI: %.4f%%\n', error_ellip_EI);
fprintf('Error in GJ: %.4f%%\n', error_ellip_GJ);

%% Function to optimize the rectangular cross-section
function [b, h, t, error_EI, error_GJ] = optimize_rectangular_cross_section(EI_target, GJ_target, E, G, lb_rec, ub_rec, supports);
    % Error function to minimize
    errorFunc_rect = @(x) abs(log(calculate_EI_rect(x(1), x(2), x(3), E, supports) / EI_target)) + ...
                     abs(log(calculate_GJ_rect(x(1), x(2), x(3), G, supports) / GJ_target));

    
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
    EI_actual = calculate_EI_rect(b, h, t, E, supports);
    GJ_actual = calculate_GJ_rect(b, h, t, G, supports);

    % Calculate the percentage errors in EI and GJ
    error_EI = ((EI_actual - EI_target) / EI_target) * 100;
    error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
end

%% Function to optimize the circular cross-section
function [r, t, error_EI, error_GJ] = optimize_circular_cross_section(EI_target, GJ_target, E, G, lb_cir, ub_cir, supports);
    errorFunc_circ = @(x) (calculate_EI_circ(x(1), x(2), E, supports) - EI_target)^2 + ...
                           (calculate_GJ_circ(x(1), x(2), G, supports) - GJ_target)^2;
    
    initial_guess_circ = [0.05, 0.001]; % Initial guesses: r, t

    % Bounds: ensure r > t, and practical limits
    lb = lb_cir* 0.0254;  % Smallest reasonable values
    ub = ub_cir * 0.0254;   % r ≤ 3 inches, t ≤ 0.01 meters

    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
    optimal_params_circ = fmincon(errorFunc_circ, initial_guess_circ, [], [], [], [], lb, ub, [], options);
    
    r = optimal_params_circ(1);
    t = optimal_params_circ(2);

    EI_actual = calculate_EI_circ(r, t, E, supports);
    GJ_actual = calculate_GJ_circ(r, t, G, supports);
    
    error_EI = ((EI_actual - EI_target) / EI_target) * 100;
    error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
end

%% Function to optimize the half-elliptical cross-section
function [a, b, t, error_EI, error_GJ] = optimize_half_elliptical_cross_section(EI_target, GJ_target, E, G, lb_eli, ub_eli, supports);
    errorFunc_ellip = @(x) (calculate_EI_ellip(x(1), x(2), x(3), E, supports) - EI_target)^2 + ...
                           (calculate_GJ_ellip(x(1), x(2), x(3), G, supports) - GJ_target)^2;
    
    initial_guess_ellip = [0.05, 0.025, 0.001]; % Initial guess: a, b, t

    % Bounds: a (major axis), b (minor axis) ≤ 3 inches, thickness reasonable
    lb = lb_eli * 0.0254;  % Smallest values
    ub = ub_eli * 0.0254;  % largest values

    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
    optimal_params_ellip = fmincon(errorFunc_ellip, initial_guess_ellip, [], [], [], [], lb, ub, [], options);
    
    a = optimal_params_ellip(1);
    b = optimal_params_ellip(2);
    t = optimal_params_ellip(3);

    EI_actual = calculate_EI_ellip(a, b, t, E, supports);
    GJ_actual = calculate_GJ_ellip(a, b, t, G, supports);
    
    error_EI = ((EI_actual - EI_target) / EI_target) * 100;
    error_GJ = ((GJ_actual - GJ_target) / GJ_target) * 100;
end


%% Function to calculate EI for a rectangular section
function EI = calculate_EI_rect(b, h, t, E, supports)
    b_inner = b - 2 * t; 
    h_inner = h - 2 * t;
    if b_inner <= 0 || h_inner <= 0
        EI = Inf; % Penalize infeasible values
    else
        EI = (E / 12) * (b * h^3 - b_inner * h_inner^3)*supports;
    end
end

%% Function to calculate GJ for a rectangular section
function GJ = calculate_GJ_rect(b, h, t, G, supports)
    b_inner = b - 2 * t;  
    h_inner = h - 2 * t;
    if b_inner <= 0 || h_inner <= 0
        GJ = Inf; % Penalize infeasible values
    else
        GJ = G * ( (1/12*b*h*(b^2+h^2)) - (1/12*b_inner*h_inner*(b_inner^2+h_inner^2)) )*supports;
    end
end

%% Function to calculate EI for a circular section
function EI = calculate_EI_circ(r, t, E, supports)
    r_inner = r - t;
    if r_inner <= 0
        EI = Inf; % Penalize infeasible values
    else
        EI = (E * pi / 4) * (r^4 - r_inner^4)*supports;
    end
end

%% Function to calculate GJ for a circular section
function GJ = calculate_GJ_circ(r, t, G, supports)
    r_inner = r - t;
    if r_inner <= 0
        GJ = Inf; % Penalize infeasible values
    else
        GJ = (G * pi / 2) * (r^4 - r_inner^4)*supports;
    end
end


%% Function to calculate EI for a half-elliptical section
function EI = calculate_EI_ellip(a, b, t, E, supports)
    a_inner = a - t;
    b_inner = b - t;
    if a_inner <= 0 || b_inner <= 0
        EI = Inf;
    else
        % Second moment of area for a half-ellipse about the flat base
        I_outer = (pi / 8 - 8/(9*pi)) * (a * b^3);
        I_inner = (pi / 8 - 8/(9*pi)) * (a_inner * b_inner^3);
        EI = E * (I_outer - I_inner)*supports;
    end
end

%% Function to calculate GJ for a half-elliptical section
function GJ = calculate_GJ_ellip(a, b, t, G, supports)
    a_inner = a - t;
    b_inner = b - t;
    if a_inner <= 0 || b_inner <= 0
        GJ = Inf;
    else
        % Approximate polar moment of inertia for a half-ellipse
        J_outer = (pi / 8 - 8/(9*pi)) * (a * b^3 + a^3 * b);
        J_inner = (pi / 8 - 8/(9*pi)) * (a_inner * b_inner^3 + a_inner^3 * b_inner);
        GJ = G * (J_outer - J_inner)*supports;
    end
end