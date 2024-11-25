clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';


syms R k0 Ea T B
%Values
R_val = 8.314; %J*K^-1 * mol^-1
Ea_val = 40 * 10^3; %J/mol
T_val = 289; %K
B_val = 0.8; 
k0_val = 1.5 * 10^6; %s^-1

%Uncertatinty
dEa_val = 10 * 10^3; %J/mol
dT_val = 5; %K
dB_val = 0.1; 
dk0_val = 0.5 * 10^6; %s^-1

k = k0*exp(-(Ea/(R*T))^B);

%Partial Derivatives
dk_dk0 = diff(k, k0);
dk_dEa = diff(k, Ea);
dk_dT = diff(k, T);
dk_DB = diff(k, B);

%Plug values into partials
dk_dk0_val= (subs(dk_dk0, {R, k0, Ea, T, B}, {R_val, k0_val, Ea_val, T_val, B_val}));
dk_dEa_val = (subs(dk_dEa, {R, k0, Ea, T, B}, {R_val, k0_val, Ea_val, T_val, B_val}));
dk_dT_val = (subs(dk_dT, {R, k0, Ea, T, B}, {R_val, k0_val, Ea_val, T_val, B_val}));
dk_dB_val = (subs(dk_DB, {R, k0, Ea, T, B}, {R_val, k0_val, Ea_val, T_val, B_val}));

%Nominal Value
k_nom = k0_val*exp(-(Ea_val/(R_val*T_val))^B_val);

%Statistical Uncertainty
dk = sqrt((dk_dk0_val * dk0_val)^2 + (dk_dEa_val * dEa_val)^2 + (dk_dT_val * dT_val)^2 + (dk_dB_val * dB_val)^2);

%Results
fprintf('The Nominal Value is %.3e\n', k_nom);
fprintf('The Statistical Uncertainty is %.3e\n', dk);

