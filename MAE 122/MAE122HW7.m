clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%% Problem 1
T = 16; %sec
omega = 2*pi/T;
g = 9.81;
gamma = 0.78;
hdeep = [1 4]; %m


%Find H from Hdeep and T
%Equation H(i) = Hdeep * (cgdeep / (2*cg(i)))^0.5 is rearranged using
%shallow approximation for cg and H = 0.78h.

for i = 1:2
   f = @(h) gamma*h - hdeep(i) * sqrt( (g*T/(4*pi)) / sqrt(g*h) );
   h_guess = 10;  % Initial guess
   h(i) = fzero(f, h_guess);
   H(i) = h(i) * 0.78;
end

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
u(i) = omega * H(i) / 2 * (cosh(ksoln(i) * h(i)) / sinh(ksoln(i) * h(i)));  % Horizontal velocity near surface (z=0)
uorbital(i) = omega * (H(i) / 2); %Amplitude of u velocity equation assuming shallow water condition (tanh(x) --> x) 

end
end

  % Display results
    for i = 1:2
    fprintf('For the case of T = 16 sec and Hdeep = %.4f m\n', hdeep(i));
    fprintf('Wave height at breaking = %.4f m\n', H(i));
    fprintf('Phase speed, c, at breaking location = %.2f m\n', c(i));
    fprintf('Maximum horizontal water velocity at breaking location %.2f m\n', max(u(i)));
    fprintf('Horizontal orbital amplitude at breaking location = %.4f m\n', uorbital(i));
    end


