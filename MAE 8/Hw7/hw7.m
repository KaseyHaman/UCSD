clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 7;

%%Problem 1
load('Terrain.mat');

figure(1);
surf(x, y, altitude'); shading interp;
xlabel('x'); ylabel('y'); zlabel('altitude(x,y)');
title('Altitude demo');
counter_lmax = 0;
snow_cover = 0;
hold on;
for n = 2:length(x)-1
     for m = 2:length(y)-1
%          nb(1) = altitude(n-1, m+1);
%          nb(2) = altitude(n-1, m);
%          nb(3) = altitude(n-1, m-1);
%          nb(4) = altitude(n, m-1);
%          nb(5) = altitude(n+1, m-1);
%          nb(6) = altitude(n+1, m);
%          nb(7) = altitude(n+1, m+1);
%          nb(8) = altitude(n+1, m+1);
          nb = altitude(n-1:n+1, m-1:m+1);
          nb = nb(:);
         if altitude(n, m) >= max(nb)
             counter_lmax = counter_lmax + 1;
             x_localmax(counter_lmax) = x(n);
             y_localmax(counter_lmax) = y(m);
             altitude_localmax(counter_lmax) = altitude(n, m);
         elseif altitude(n, m) >= 1100
             snow_cover = snow_cover + 1;
             x_snow(snow_cover) = x(n);
             y_snow(snow_cover) = y(m);
             altitude_snow(snow_cover) = altitude(n, m);
         end
     end
end
plot3(x_localmax, y_localmax, altitude_localmax, 'ro', 'Markersize', 10, 'MarkerFacecolor', 'r');
plot3(x_snow, y_snow, altitude_snow, 'go', 'Markersize', 4, 'MarkerFacecolor', 'g');
legend('Terrain', 'Local Peaks', 'Snow Cover', 'location', 'best')

p1a = counter_lmax; %Double check
p1b = x_localmax;
p1c = y_localmax;
p1d = altitude_localmax;
p1e = x_snow;
p1f = y_snow;
p1g = altitude_snow;
p1h = 'See Figure 1';

%%Problem 2

load('matB.mat');
Bsum = 0;
p2a = 0;

sum1 = 0;
for n = 1:length(matB)
   for m = 1:length(matB)
       if n<m
           sum1 = sum1 + matB(n, m);
       elseif n>=m
           continue;
       end
   end
end
p2a = sum1;

prod2 = 1;
for n = 1:length(matB)
   for m = 1:length(matB)
       if n>m
           prod2 = prod2 * matB(n, m);
       elseif n<=m
           continue;
       end
   end
end
p2b = prod2;

sum3 = 0;
for n = 1:length(matB)
   for m = 1:length(matB)
       if 2*n == m
           continue;
       else
           sum3 = sum3 + matB(n, m);
       end
   end
end
p2c = sum3;

prod4 = 1;
for n = 1:length(matB)
   for m = 1:length(matB)
       if matB(n, m) > n
           continue;
       else
           prod4 = prod4 * matB(n, m);
       end
   end
end
p2d = prod4;


%%Problem 3
p3a = evalc('help car');
Tf = 60;
dt = 10;
[T, X, U] = car(Tf, dt);
p3b = T;
p3c = X;
p3d = U;
Tf = 60;
dt = 1;
[A, B, C] = car(Tf, dt);
p3e = A;
p3f = B;
p3g = C;

figure(2);
subplot(2, 1, 1); hold on;
plot(T, X, '-r');
plot(A, B, '-bo');
ylabel('Distance (m)');
xlabel('Time (s)');
title('Motion of the Car');
hold off;
subplot(2, 1, 2); hold on;
plot(T, U, '-r');
plot(A, C, '-bo');
ylabel('Velocity (m/s)');
xlabel('Time (s)');
title('Motion of the Car');
hold off;
p3h = 'See figure 2';

%%Problem 4
p4a = evalc('help rocket');
p4b = evalc('help gravity');
p4c = evalc('help thrust');
p4d = evalc('help rocket>gravity');
p4e = evalc('help rocket>thrust');
Tf = 120; dt = 0.1;
[T, Z, W] = rocket(Tf, dt);
p4f = Z(end);
p4g = W(end);
figure(3);
subplot(2, 1, 1); 
plot(T, Z, '-r');
ylabel('Distance (m)');
xlabel('Time (s)');
title('Change in Altitude of Rocket Over Time');
subplot(2, 1, 2); 
plot(T, W, '-b');
ylabel('Velocity (m/s)');
xlabel('Time (s)');
title('Change in Velocity of Rocket Over Time');
p3h = 'See figure 2';

