clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1

%Part A
TP = readmatrix('TP.dat');
x = TP(:, 1);
y = TP(:, 2);

scatter(x,y, 10,'r');
grid on;
axis on;
xlabel('Temperature')
ylabel('Power')

%Part B

%i) y = Ax + a --> Ax = b ---> x = A^-1 *b ........NOTE: the x values in
%this case are the c constants (what we're solving for)
xi = [x ones(size(x, 1), 1)];
linearcoeff = xi \ y;

%ii) y = Ax^2 + bx + a --> Ax = b ---> x = A^-1 *b
xii = [x.^2 x ones(size(x, 1), 1)];
quadraticcoeff = xii \ y;

%iii) y = Ax^3 + bx^2 + cx + a --> Ax = b ---> x = A^-1 *b
xiii = [x.^3 x.^2 x ones(size(x, 1), 1)];
cubiccoeff = xiii \ y;

%iv) y = Ax^4 + bx^3 + cx^2 + dx + a --> Ax = b ---> x = A^-1 *b
xiv = [x.^4 x.^3 x.^2 x ones(size(x, 1), 1)];
quarticcoeff = xiv \ y;

%Part C
hold on;
plot(x, xi*linearcoeff, '-b', 'LineWidth', 1.5)
plot(x, xii*quadraticcoeff, '-k', 'LineWidth', 1.5)
plot(x, xiii*cubiccoeff, '-g', 'LineWidth', 1.5)
plot(x, xiv*quarticcoeff, '-m', 'LineWidth', 1.5)
title('My fits')

%Part D
linearRMSE = sqrt(mean((y-xi*linearcoeff).^2));
quadraticRMSE = sqrt(mean((y-xii*quadraticcoeff).^2));
cubicRMSE = sqrt(mean((y-xiii*cubiccoeff).^2));
quarticRMSE = sqrt(mean((y-xiv*quarticcoeff).^2));

if linearRMSE > quadraticRMSE
    a1 = quadraticRMSE;
    a1p = 'quadraticRMSE';
else
    a1 = linearRMSE;
    a1p = 'linearRMSE';
end

if cubicRMSE > quarticRMSE
    a2 = quarticRMSE;
    a2p = 'quarticRMSE';
else
    a2 = cubicRMSE;
    a2p = 'cubicRMSE';
end

if a1 > a2
    a = a2p;
    a0 = a2;
else
    a = a1p;
    a0 = a1;
end

sprintf('The lowest RMSE value belonged to %s at %d', a, a0)

figure(2);
ylinear = polyfit(x, y, 1);
yquadratic = polyfit(x, y, 2);
ycubic = polyfit(x, y, 3);
yquartic = polyfit(x, y, 4);
scatter(x,y, 10,'r');
grid on;
axis on;
xlabel('Temperature')
ylabel('Power')
hold on;
plot(x, xi*ylinear', '-b', 'LineWidth', 1.5)
plot(x, xii*yquadratic', '-k', 'LineWidth', 1.5)
plot(x, xiii*ycubic', '-g', 'LineWidth', 1.5)
plot(x, xiv*yquartic', '-m', 'LineWidth', 1.5)
title('MATLABs polyfit')






%%Problem 2
clear all;

%Part A
tV = readmatrix('tV.dat');
x = tV(:, 1);
y = tV(:, 2);

figure(3);
scatter(x,y, 10,'r');
grid on;
axis on;
xlabel('Time')
ylabel('Voltage')
title('Voltage and Time fit')

%Part B

%On paper work.

%Part C

yi = [sum(x.*y); sum(y.*sin(2*pi.*x))];
xi = [sum(x.^2) sum(x.*sin(2*pi.*x)); sum(x.*sin(2*pi.*x)) sum(sin(2*pi.*x).^2)];
sincoeff = xi \ yi;

% %Alternatively) y = ax + b(sin(2*pi*x) --> Ax = b ---> x = A^-1 *b
% % Note: This uses MATLAB's built in least squares solver.
% xi = [x sin(2*pi*x)];
% sincoeff = xi \ y;

sprintf('The sin coefficients are %d and %d', sincoeff(1), sincoeff(2))

%Part D
hold on;
xcoeff = [x sin(2*pi*x)];
plot(x, xcoeff*sincoeff, '-k', 'LineWidth', 1.5)

%Part E

% Define the model function: V(t) = c1 * x + c2 * sin(2*pi*x)
func = @(c, x) c(1) * x + c(2) * sin(2 * pi * x);

% Coefficient Guess
coeffguess = [1, 1];

% Perform the curve fitting using lsqcurvefit
options = optimoptions('lsqcurvefit', 'Display', 'off');
lsqcurvefitcoeff = lsqcurvefit(func, coeffguess, x, y, [], [], options);

%Plot the lsqcurvefit curve
plot(x, xcoeff*lsqcurvefitcoeff', '--g', 'LineWidth', 1)
legend('Raw data', 'My fit', 'Matlabs fit')





%%Problem 3
clear all;

%Part A
xy = readmatrix('xy.dat');
x = xy(:, 1);
y = xy(:, 2);

figure(4);
scatter(x,y, 10,'b');
grid on;
axis on;
xlabel('x')
ylabel('y')
title('Exponential Regression Model')

%Part B

%Completed on paper

%Part C

%Note that these linear coefficients are the solution to y' = b'x + a' 
%where y' = ln(y - 273), b' = b*ln(10) and a' = ln(a) from the equation
% y(x) = a * 10^(b*x) + 273


% Linearizing: y(x) = a * 10^(b*x) + 273
% We need to calculate ln(y - 273)
y_shifted = y - 273; % Shift y
ylinearized = log(y_shifted); % Take the natural log of (y - 273)

xi = [x ones(size(x, 1), 1)];
loglinearcoeff = xi \ ylinearized;

%Solve for coefficient b
linearcoeff(1) = loglinearcoeff(1) / log(10);

%Solve for coefficient a
linearcoeff(2) = exp(loglinearcoeff(2));

%Plot the best fit approximation

hold on;
plot(x, linearcoeff(2).*10.^(linearcoeff(1).*x)+273, '-k', 'LineWidth', 1.5)
