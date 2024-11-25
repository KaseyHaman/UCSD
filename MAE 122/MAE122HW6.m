clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1
T = 16; %sec
w = 2*pi/T;
g = 9.81;
h = [10 20 50 100 400 800];
soln = zeros(size(h));

for i = 1:6
    %Calculate k given dispersion relation (w^2 - gktanh(kh) = 0)
    %Bisections Method to find soln;
func = @(k) w^2-g*k*tanh(k*h(i));
a = 0; 
b = 0.10; %lower and upper bounds
if func(a) * func(b) > 0
    fprintf('Error with bounds at h=%d', h(i))
    break
end
while (b-a)/2 > 1e-7
c = (a+b)/2;
if func(a)*func(c) < 0
    b = c;
else
    a = c;
end
    ksoln(i) = (b+a)/2;

end
end

clear c;
%Calculate phase and group speeds
for i = 1:6
    c(i) = w/ksoln(i);
    cg(i) = 0.5 * (w / ksoln(i)) * (1 + (2 * ksoln(i) * h(i)) / sinh(2 * ksoln(i) * h(i)));
    lamda(i) = 2*pi / ksoln(i);
end

%Calculate Deep and Shallow Water Limits
for i = 1:6
kdeep = w^2/g;
cdeep(i) = sqrt(g/kdeep);
kshallow = sqrt(w^2/(g*h(i)));
cshallow(i) = sqrt(g*h(i));

cgroupshallow(i) = sqrt(g*h(i));
cgroupdeep(i) = c(i)/2;
end

figure;
subplot(3, 1, 1)
plot (h, lamda, 'b-o')
hold on;
title('Wavelength Over Varying Depth')
ylabel('Wavelength (m)');
xlabel('Depth (m)');

subplot(3, 1, 2)
plot (h, c, 'k-o', 'LineWidth', 1.2)
hold on;
title('Speed Over Varying Depth')
ylabel('Speed (m/s)');
xlabel('Depth (m)');
plot(h, cdeep, 'g--')
plot(h, cshallow, 'm--')
legend('Phase Speed',  ...
    'Deep Approx.', 'Shallow Approx.',  'Location', 'Northeast');

subplot(3, 1, 3)
plot(h, cg, 'k-o', 'LineWidth', 1.2)
hold on;
title('Group Speed Over Varying Depth')
ylabel('Speed (m/s)');
xlabel('Depth (m)');
plot(h, cgroupdeep, 'g--')
plot(h, cgroupshallow, 'm--')
legend('Group Speed', ...
    'Deep Approx.', 'Shallow Approx.',  'Location', 'Northeast');


%% Problem 2
clear all;
h = [10, 30];      % m
e = 0.6;           % efficiency constant
T = 8;             % s
Hdeep = 3;             % m
p = 1029;          % kg/m^3 seawater
desiredP = 500;    % watts
g = 9.81;          % m/s^2
omega = 2 * pi / T;   %rad/s

%Solve for deep water group velocity and wavelength
cgdeep = 1/2*(g*T/(2*pi));
lamdadeep = cgdeep*T*2;

%Calculate k values at different h values
for i = 1:2
    %Calculate k given dispersion relation (w^2 - gktanh(kh) = 0)
    %Bisections Method to find soln;
func = @(k) omega^2-g*k*tanh(k*h(i));
a = 0; 
b = 10; %lower and upper bounds
if func(a) * func(b) > 0
    fprintf('Error with bounds at h=%d', h(i))
    break
end
while (b-a)/2 > 1e-7
d = (a+b)/2;
if func(a)*func(d) < 0
    b = d;
else
    a = d;
end
    ksoln(i) = (b+a)/2;

%Solve for other constant(s)
lamda(i) = 2*pi/ksoln(i); %wavelength
c(i) = omega / ksoln(i);
cg(i) = c(i)/2 * ( 1 + (2*ksoln(i)*h(i)) / sinh(2*ksoln(i)*h(i)) );

%Solve H with equation from lecture
H(i) = Hdeep * (cgdeep / (2*cg(i)))^0.5;

u(i) = omega * H(i) / 2 * (1 / sinh(ksoln(i) * h(i)));  % Horizontal velocity near sea bed (z=-h)
w(i) = omega * H(i) / 2;  % Vertical velocity near surface (z=0)

% Diameter Calculations
Dhorizontal(i) = 2 * sqrt(desiredP / (e * 2 * pi * 0.5 * p * (u(i)^3)));  % Horizontal cylinder diameter
Dvertical(i) = 2 * sqrt(desiredP / (e * 2 * pi * 0.5 * p * (w(i)^3)));  % Vertical cylinder diameter

end
end


    % Display results
    for i = 1:2
    fprintf('For depth = %d m:\n', h(i));
    fprintf('  Wavelength = %.2f m\n', lamda(i));
    fprintf('  Required diameter for vertical cylinder: %.2f m\n', Dvertical(i));
    fprintf('  Required diameter for horizontal cylinder: %.2f m\n\n', Dhorizontal(i));
    end
