function [T,X,Y,Z,U,V,W] = soccer(omega)
% Soccer computes the trajectory of a soccer ball based on the spinning rate omega in 3D
% Cartesian coordiate using Euler method. Inputs are the drag coefficient. 
% Output vectors are time (T), position (X,Y,Z) and velocity components (U,V,W).
% Call format: [T,X,Y,Z,U,V,W] = soccer(omega)


% Set up parameters
r = 0.11; A = pi*r^2; m = 0.4; rho = 1.2; 
g = 9.81; dt = 1/1000; Cm = -0.6;

% Initialize the kick
n = 1;
T(n) = 0;
X(n) = 0;
Y(n) = 0;
Z(n) = 0;
U(n) = 10;
V(n) = -20;
W(n) = 10;
 
% Advance the governing equation via Euler method
while Z(n) >= 0 
    coeff = Cm*rho*A*r/2/m;
    U(n+1) = U(n) + coeff*omega*V(n)*dt;
    V(n+1) = V(n) - coeff*omega*U(n)*dt;
    W(n+1) = W(n) - g*dt;
    X(n+1) = X(n) + U(n)*dt;
    Y(n+1) = Y(n) + V(n)*dt;
    Z(n+1) = Z(n) + W(n)*dt;
    T(n+1) = T(n) + dt;         
    n = n+1;
end

% Remove the last element where the soccer is below ground
X(end) = []; Y(end) = []; Z(end) = []; 
U(end) = []; V(end) = []; W(end) = [];
T(end) = [];
end %function soccer