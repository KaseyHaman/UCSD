clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1

%Condition 1

%initialize x y and z.
x(1) = 3; y(1) = 2; z(1) = 1;

%initialize functions
f0 = @(x, y, z) 1/4*(sin(x + y + atan(z/2)));
f1 = @(x, y, z) 1/4*(cos(x + y + atan(z/2)));
f2 = @(x, y, z) 1 + 1/4*(cos(x + atan(z/2)));

n = 100; %iterations
i = 1;
while i <= n
xstep = f0(x(i), y(i), z(i));
ystep = f1(x(i), y(i), z(i));
zstep = f2(x(i), y(i), z(i));
i = i+1;
x(i) = xstep;
y(i) = ystep;
z(i) = zstep;
end

%Next 3 iterations
condition1_x1_iterations = x(2:4)
condition1_x2_iterations = y(2:4)
condition1_x3_iterations = z(2:4)

%Estimate of solution
condition1_x1_estimate = x(end)
condition1_x2_estimate = y(end)
condition1_x3_estimate = z(end)

%Condition 2

%initialize x y and z.
x(1) = 300; y(1) = 200; z(1) = 100;

%initialize functions
f0 = @(x, y, z) 1/4*(sin(x + y + atan(z/2)));
f1 = @(x, y, z) 1/4*(cos(x + y + atan(z/2)));
f2 = @(x, y, z) 1 + 1/4*(cos(x + atan(z/2)));

n = 100; %iterations
i = 1;
while i <= n
xstep = f0(x(i), y(i), z(i));
ystep = f1(x(i), y(i), z(i));
zstep = f2(x(i), y(i), z(i));
i = i+1;
x(i) = xstep;
y(i) = ystep;
z(i) = zstep;
end

%Next 3 iterations
condition2_x1_iterations = x(2:4)
condition2_x2_iterations = y(2:4)
condition2_x3_iterations = z(2:4)

%Estimate of solution
condition2_x1_estimate = x(end)
condition2_x2_estimate = y(end)
condition2_x3_estimate = z(end)


%%Problem 2
clear x;
t(1) = 1; t(2) = 3; t(3) = 4; t(4) = 5; t(5) = 7;
x(1) = -0.9; x(2) = 1.7; x(3) = 2.1; x(4) = 0.2; x(5) = -0.7;

%Initialize parts of matrix (r = row, c = column)
r1c1 = sum(t.^4); r1c2 = sum(t.^3); r1c3 = sum(t.^2); b1 = sum(t.^2.*x);
r2c1 = sum(t.^3); r2c2 = sum(t.^2); r2c3 = sum(t); b2 = sum(t.*x);
r3c1 = sum(t.^2); r3c2 = sum(t); r3c3 = size(t, 2); b3 = sum(x);
%Initalize matrix
A = [r1c1 r1c2 r1c3; r2c1 r2c2 r2c3; r3c1 r3c2 r3c3]; b = [b1; b2; b3];
%Because Ax = b, x = A\b
soln = A\b;
a = soln(1);
b = soln(2);
c = soln(3);
%Write function
func = @(x) a.*x.^2 + b.*x + c;
%Plot points
tdomain = 0:0.5:10;
figure(1), hold on;
% cs = 'krbgmckrbgm';
 plot(tdomain,func(tdomain), 'k','LineWidth',1);
 plot(t,x, 'bo');
%  legend('Function Evalutation at x', 'Location', 'best')
xlabel('t'); ylabel('function at t');
title('Least-squares Quadratic Function');
box on; grid on;
set(gca,'FontSize',10)


