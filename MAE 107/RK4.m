function [yk_values,t_values] = RK4(f, y, z, h)
%This is a Runge Kutta 4 function designed to find the root of a given
%numerical function.
%Note that this function includes the fixed point method used for solving
%the root of a function given 3 variables with 1 unknown.

%Run Runge Kuta 4 loop
 i = 0; y(1) = y; z(1) = z;
    for t = 0:h:6
        i = i + 1;
        k1 = h*f(z, t, y(i));
        k2 = h*f(z, t+h/2, y(i)+k1/2);
        k3 = h*f(z, t+h/2, y(i)+k2/2);
        k4 = h*f(z, t+h, y(i)+k3);
        y(i+1) = y(i) + 1/6*(k1+2*(k2+k3)+k4);
        z = f(z, t, y(i));
    end
yk_values = y(1:end-1);
t_values = 0:h:6;


end