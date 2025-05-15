function [mu] = sutherland(T)
%Function to compute dynamic viscosity from empirical data using sutherland
%equation

%Define constants
mu_o = 1.735e-5; %N*s/m^2
T_o = 288; %K
S_1 = 110.4; %K

mu = mu_o*(T/T_o).^(3/2).*((T_o + S_1)./(T + S_1));

end

