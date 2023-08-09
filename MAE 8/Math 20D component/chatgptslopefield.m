function chatgptslopefield(ODE, xrange, yrange, stepsize)
% This function generates a slope field of a given ODE
% 
% Inputs:
%   ODE: a function handle representing the ODE, in the form y' = f(x,y)
%   xrange: a 1x2 array representing the range of x values to plot
%   yrange: a 1x2 array representing the range of y values to plot
%   stepsize: a scalar representing the stepsize for the slope field
%
% Output:
%   A slope field plot of the ODE

% Define the x and y values for the slope field
[x, y] = meshgrid(xrange(1):stepsize:xrange(2), yrange(1):stepsize:yrange(2));

% Calculate the slope at each point using the ODE function
dy = ODE(x, y);
dx = ones(size(dy));

% Normalize the slope to the stepsize
slope = stepsize*dy./dx;

% Plot the slope field
quiver(x, y, dx, dy, 'k')
axis tight equal;
xlabel('x');
ylabel('y');
end