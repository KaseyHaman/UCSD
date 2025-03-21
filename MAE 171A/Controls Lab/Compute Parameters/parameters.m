% PARAMETERS - Belongs to MAELAB.M for MAE171a
% MatLab template file for specifying 2DOF model parameters

% Written by R.A. de Callafon, Dept. of MAE, UCSD (2001-2021)
% Report errors in this software to <callafon@ucsd.edu>

% Define file paths
files_Mass1Still = {'Mass1Still1.mat', 'Mass1Still2.mat', 'Mass1Still3.mat', 'Mass1Still4.mat', 'Mass1Still5.mat'};
files_Mass2Still = {'Mass2Still1.mat', 'Mass2Still2.mat', 'Mass2Still3.mat', 'Mass2Still4.mat', 'Mass2Still5.mat'};
files_TwoDOF = {'TwoDOF1.mat', 'TwoDOF2.mat', 'TwoDOF3.mat', 'TwoDOF4.mat', 'TwoDOF5.mat'};

F = 0.5;  % Example force
n = 1;    % Example number

% Call the compute_parameters function to get the average values
[AvgK1, AvgK2, AvgM1, AvgM2, AvgD1, AvgD2] = compute_parameters(files_Mass1Still, files_Mass2Still, files_TwoDOF, F, n);

% Use the computed parameters for the 2DOF model
test = 2;  % This variable no longer needs to index the averaged results
Mass1 = AvgM1;
Mass2 = AvgM2;
Damp1 = AvgD1;
Damp2 = AvgD2;
Spr1 = AvgK1;
Spr2 = AvgK2;
Spr3 = Spr1 + Spr2;
% Mass1 = 1.61e-06;
% Mass2 = 1.43e-06;
% Damp1 = 6.48e-06;
% Damp2 = 3.20e-06;
% Spr1 = 4.91e-04;
% Spr2 = 3.56e-04;
% Spr3 = Spr1 + Spr2;

% Since AvgM1, AvgM2, AvgD1, AvgD2, AvgK1, and AvgK2 are scalars, you can assign them directly
m1 = Mass1;         % Mass m1 (scalar)
d1 = Damp1;         % Damping that connects m1 to ground (scalar)  
k1 = Spr1;          % Spring that connects m1 to ground (scalar)

m2 = Mass2;         % Mass m2 (scalar)
d2 = Damp2;         % Damping that connects m2 to ground (scalar)
k2 = Spr2;          % Spring that connects m1 and m2 (scalar)

% Now you can use the values of m1, m2, d1, d2, k1, k2 for your 2DOF model
