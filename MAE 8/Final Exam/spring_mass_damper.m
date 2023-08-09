function [T, X, V] = spring_mass_damper(c)
%This function takes the frictional damper factor, c, and outputs
%the time, displacement and velocity of mass for a 2 second interval.


% Set up parameters
m = 0.1; k = 40; dt = 0.001; 

% Initialize
n = 1;
T(n) = 0;
X(n) = 0.1;
V(n) = 10;

while T(n) <= 2
    V(n+1) = V(n) - (k/m*X(n) + c/m*V(n))*dt;
    X(n+1) = X(n) + V(n+1)*dt;
    T(n+1) = T(n) + dt;
    n = n+1;
end
V(end) = [];
X(end) = [];
T(end) = [];

end

