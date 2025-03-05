clear; close all; clc;
%% Dogbone Tensile Test Calculations
% Given Material Properties
sigma_ult = 65e6;   % Ultimate tensile strength (Pa)
Emin = 0.0420e9;         % Elastic modulus (Pa)
Emax = 3.3e9;
nu = 0.35;         % Poisson's ratio

% Dogbone Specimen Dimensions
w = 13e-3;         % Width of narrow section (m)
t = 3.2e-3;        % Thickness (m)
L_g = 50e-3;       % Gage length (m)

% Testing Equipment Specifications
load_cell_capacity = 30e3;  % Load cell max force (N)
load_frame_displacement = 100; % Max displacement (mm)

% Camera Specifications
resolution = 2500; % iPhone estimated horizontal resolution (pixels)
fov = 250;        % Field of view in mm

%% Problem 1: Maximum Expected Load
A = w * t; % Cross-sectional area (m^2)
F_max = sigma_ult * A; % Maximum expected load (N)

fprintf('Problem 1: Maximum Expected Load\n');
fprintf('Maximum Load: %.2f kN\n', F_max / 1e3);
if F_max < load_cell_capacity
    fprintf('Load cell capacity is sufficient.\n\n');
else
    fprintf('WARNING: Load cell capacity may be exceeded!\n\n');
end

%% Problem 2: Maximum Expected Displacement
strain_max = sigma_ult / Emin; % Maximum strain (unitless)
strain_min = sigma_ult / Emax; % Maximum strain (unitless)
delta_max = strain_max * L_g; % Maximum displacement (m)
delta_min = strain_min * L_g; % Maximum displacement (m)

delta_max_mm = delta_max * 1e3; % Convert to mm
delta_min_mm = delta_min * 1e3; % Convert to mm
fprintf('Problem 2: Maximum Expected Displacement\n');
fprintf('Maximum Displacement: %.2f mm\n', delta_max_mm);
fprintf('Minimmum Displacement: %.2f mm\n', delta_min_mm);
if delta_max_mm < load_frame_displacement
    fprintf('Load frame displacement is sufficient.\n\n');
else
    fprintf('WARNING: Load frame displacement may be exceeded!\n\n');
end

%% Problem 3: Strain and Displacement Resolution using DIC
pixel_size = fov / resolution; % mm per pixel
strain_resolution = pixel_size / (L_g*10^3); % Strain resolution (unitless)

fprintf('Problem 3: Strain & Displacement Resolution\n');
fprintf('Pixel Size: %.4f mm/pixel\n', pixel_size);
fprintf('Strain Resolution: %.5f (or %.2f%%)\n', strain_resolution, strain_resolution * 100);
