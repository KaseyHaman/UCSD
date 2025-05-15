function [mu] = sutherland(T)
%Solves Sutherland Equation

%Define Variables for air 
mu0 = 1.735*10^-5; %N*s / m^2
S1 = 110.4; %k
T0 = 288; %k

%Solve Equation
mu = mu0.*(T/T0).^(3/2).*(T0 + S1)./(T+S1);

end