close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';
hw_num = 'final_practice';
form = 'A';

% %%Problem 1
% theta = -10:0.1:350;
% x = cosd(4.*theta).*cosd(theta);
% y = cosd(4.*theta).*sind(theta);
% z = cosd(8.*theta);
% figure (1); hold on; view(45, 60); grid on; box on;
% plot3(x, y, z, '-m', 'linewidth', 4)
% plot3(x(1:101:3601), y(1:101:3601), z(1:101:3601), 'go', 'markerfacecolor', 'g')


%%Problem 2
rng(int8(form));
for n = 1:3
for m = 1:3
cellA{n,m}= randi(100,n,m);
structA(n).field1 = char(80+[1:(3*n)]);
structA(n).field2 = randi(100,n,m);
end
end
p2a = cellA{3, 2}(2, 2);
p2b = mean(mean(cellA{3, 3}));

for n = 1:9
    maxval(n) = max(max(cellA{n}));
end
p2c =max(maxval);
p2d = structA(2).field1(end);
p2e = structA(3).field2(2, 2);
p2f = structA;
for n = 1:3
   p2f(n).field1(3) = p2f(n).field1(3) + 32;
end


%%Problem 3
x = -20:0.1:20;
y = cos(2.*x).*tanh(x./10);
p3a = x;
p3b = y;
p3c = (x(2) - x(1)) .* (0.5 .* (y(end) + y(1)) + sum(y(2:end-1)));

countermax = 0;
ymax = 0;
for n = 2:length(y)-1
    if y(n) > y(n+1) && y(n) > y(n-1)
        countermax = countermax + 1;
        ymax(countermax) = y(n);
    end
end


countermin = 0;
ymin = 0;
for n = 2:length(y)-1
    if y(n) < y(n+1) && y(n) < y(n-1)
        countermin = countermin + 1;
        ymin(countermin) = y(n);
    end
end

crosszero = 0;
for i = 2:length(y)-1
        if (y(i)<0 && y(i+1) > 0) || (y(i)>=0 && y(i+1)< 0) || (y(i) == 0)
        crosszero = crosszero + 1;
        sp3h(crosszero) = x(n);
        sp3i(crosszero) = y(n);
    end
end

%%Problem 4
p4a = 0;
for m = 1:100
    for n = 1:100
        p4a = p4a + 1/(2^(m*n));
    end
end

savings = 0;
a = 1;
day = 0;
while savings < 10000
savings = savings + a;
a = a*3;
day = day + 1;
end
p4b = day;

a = 1;
n = 0;
suma = 0;
while abs(pi - 4*a) > 10e-3
    suma = suma + ((-1)^n)/(2*n+1);
    a = suma;
    n = n+1;
end
p4c = n-1;
%%Problem 5

