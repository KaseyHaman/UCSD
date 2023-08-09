clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';


%%Problem 5
func = @(x) sin(pi*x) - 1/2*cos(2*pi*x);
%Left Endpoint rule
for k = 1:4
    n(k)= 10^k;
    leftendpoint(func, 0, 3, n(k));
    LE(k) = ans;
end
%Trapezoid rule
for k = 1:4
    n(k)= 10^k;
    trap(func, 0, 3, n(k));
    Trap(k) = ans;
end
%Corrected Trapezoid rule
for k = 1:4
    n(k)= 10^k;
    ctrap(func, 0, 3, n(k));
    CTrap(k) = ans;
end
%Integral of y solved by hand
yint = (-4*cos(pi*3)-sin(2*pi*3) ) / (4*pi)  -  (-4*cos(pi*0)-sin(2*pi*0) ) / (4*pi);

%Find the errors
for k = 1:4
    errorLE(k) = abs(yint - LE(k));
    errorTrap(k) = abs(yint - Trap(k));
    errorCTrap(k) = abs(yint - CTrap(k));
end
%Graph the errors
logplot(1, n, errorLE, errorTrap, errorCTrap)





%%Problem 6
func = @(x) sin(pi*(1-x^2))/(sqrt(2+x^2));

%Estimated Truth
    n= 10^5;
    ctrap(func, 0, 3, n);
    CTrap = ans;

%Left Endpoint rule
for k = 1:4
    n(k)= 10^k;
    leftendpoint(func, 0, 3, n(k));
    LE(k) = ans;
end
%Trapezoid rule
for k = 1:4
    n(k)= 10^k;
    trap(func, 0, 3, n(k));
    Trap(k) = ans;
end

%Find the errors
for k = 1:4
    errorLE(k) = abs(CTrap - LE(k));
    errorTrap(k) = abs(CTrap - Trap(k));
end

%Plot the graph
logplot(2, n, errorLE, errorTrap, 0)




%%Problem 7
func = @(x) abs(x - sqrt(2));
%Left Endpoint rule
for k = 1:4
    n(k)= 10^k;
    leftendpoint(func, 0, 3, n(k));
    LE(k) = ans;
end
%Trapezoid rule
for k = 1:4
    n(k)= 10^k;
    trap(func, 0, 3, n(k));
    Trap(k) = ans;
end
%Corrected Trapezoid rule
for k = 1:4
    n(k)= 10^k;
    ctrap(func, 0, 3, n(k));
    CTrap(k) = ans;
end
%Integral of y solved by hand
yint = 2 - (3*2.^(3/2) - 9)/2;
%Find the errors
for k = 1:4
    errorLE(k) = abs(yint - LE(k));
    errorTrap(k) = abs(yint - Trap(k));
    errorCTrap(k) = abs(yint - CTrap(k));
end
%Graph the errors
logplot(3, n, errorLE, errorTrap, errorCTrap)
