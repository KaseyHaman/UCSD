clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 3;

%%Problem 1
a1 = 1:30;
a1(1:30) = 5;
a2 = 1:29;
a2(1:29) = -2;
a3 = diag(a1);
a4 = diag(a2, 1);
a5 = diag(a2, -1);

A = a3 + a4 + a5;
b = 1:30; b(1:30) = 0; b(1) = -2; b(30) = 1; b = b';

p1 = A\b;

%%Problem 2
x = 0.0:0.1:10;
p2a = x;
f = tanh(0.5.*(x)).^4.*(exp(1).^(-(sin(x).^2)));
p2b = f;
p2c = diff(f)./diff(x);
p2d = p2c(51);

%%Problem 3
z = -5:0.1:5;
g = ((sech(z)).^2) .* ((sin(4*z)).^4);
p3a = z;
p3b = g;
p3c = (z(2) - z(1)) * (0.5 * (g(end) + g(1)) + sum(g(2:end-1)));

%%Problem 4
matA = 1:100; matA = abs(fix(100*cos(matA))); matA = reshape(matA, 10 ,10);
p4a = matA;
p4a(matA==max(matA(:, 1:10))) = -1;
p4b = matA;
p4b(matA==max(matA(:))) = -2;
p4c = isprime(matA);
p4d = find(p4c == 1);


%%Problem 5
p5a = clock;
a = p5a(1); %year
b = p5a(2); %month
c = p5a(3); %day 
d = p5a(4); %hour
e = p5a(5); %min
f = p5a(6); %second

p5b = sprintf('%4i:%02i:%02i', a, b, c);
p5c = sprintf('%02i:%02i:%07.4f', d, e, f);
p5d = sprintf('%s\b\b\b\b\b', p5c);
p5e = sprintf('%s %s', p5b, p5d);

%%Problem 6 
theta = 1:1:360;
x = 16*sind(theta).^3;
y = 13*cosd(theta)-5*cosd(2*theta)-2*cosd(3*theta)-cosd(4*theta);
figure (1);
plot(x, y, '-r', 'Linewidth', 5);
hold on;
plot(x(360), y(360), '-cd', 'Markersize', 30, 'Markerfacecolor', 'c');
grid on;
title('Two-Dimensional Parametric Curve')
xlabel('x')
ylabel('y')
set(gca, 'Fontsize', 18)
legend('Parametric Curve')
hold off;
p6a = 'See figure 1';
p6b = sum(sqrt(((diff(x).^2 + diff(y).^2))));


%%Problem 7
theta = 0:0.5:1200;
x = (1+0.25*cosd(50*theta)).*cosd(theta);
y = (1+0.25*cosd(50*theta)).*sind(theta);
z = (pi*theta)/180 + 2*sind(50 * theta);
p7a = sum(sqrt(((diff(x).^2 + diff(y).^2 + diff(z).^2))));
a = sqrt(((diff(x).^2 + diff(y).^2 + diff(z).^2)));
b = cumsum(a) ;
c = find(b <= 500);
p7b = [x(c(end)+1), y(c(end)+1), z(c(end)+1)];
figure (2); 
plot3(x, y, z, '-m', 'Linewidth', 0.5);
hold on;
plot3(x(1), y(1), z(1), '-ko', 'markersize', 10, 'Markerfacecolor', 'k');
plot3(x(end), y(end), z(end), '-r^', 'markersize', 10, 'Markerfacecolor', 'r');
plot3(x(c(end)+1), y(c(end)+1), z(c(end)+1), '-bs', 'markersize', 10, 'Markerfacecolor', 'b');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
legend('Three Dimensional Parametric Curve', 'Starting Point', ...
    'Ending Point', 'Arc Length is 500',  'Location', 'Northwest');
title('Graph of x, y and z')
hold off;
view(3);
p7c = 'See figure 2';
