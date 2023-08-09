clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 4;

%%Problem 1
n = 2:13;
A = zeros(4, 3);
A(:)= log(1./n);
A = A';
p1a = sprintf('%9.6f\n%9.6f\n%9.6f', A(:, 4));
p1b = sprintf('%7.4f %7.4f %7.4f %7.4f', A(3, :));
p1c = sprintf(['%10.7f %10.7f %10.7f %10.7f\n%10.7f %10.7f %10.7f %10.7f\n' ...
    '%10.7f %10.7f %10.7f %10.7f'], A');
p1d = sprintf(['%10.7e %10.7e %10.7e %10.7e\n%10.7e %10.7e %10.7e %10.7e\n' ...
    '%10.7e %10.7e %10.7e %10.7e'], A');

%%Problem 2
load('temperature.dat');
[r, c] = find(temperature(:, 2:13)==max(max(temperature(:, 2:13))));
p2a = c;
p2b = temperature(r, 1);

start_year = find(temperature(:, 1) == 1900);
end_year = find(temperature(:, 1) == 2000);
p2c = sum(sum(temperature(start_year:end_year, 2:13)>75));

p2d = sum(sum(temperature(start_year:end_year, 7:9)>75));

figure(1);
bar(1:12, mean(temperature(:, 2:13))); 
xlabel('Month');
ylabel('Mean Temperature(^\circF)'); 
xticks(1:12);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'});
xtickangle(45);
title('Annual Cycle of temperature in San Diego');
set(gca, 'FontSize', 16);
box on; grid on;
p2e = 'See figure 1';

p2f = 'In general, the temperature peaks in August';

output_data = [(1:12)' mean(temperature(:, 2:13))'];
save('annual_cycle.dat', 'output_data', '-ascii');
p2g = evalc('type("annual_cycle.dat")');

year = temperature(:, 1)';
annualT = mean(temperature(:, 2:13)');
figure(2); 
plot(year, annualT, '-k', 'linewidth', 2);
xlabel('Year');
ylabel('Mean Temperature(^\circF)');
title('Annual Mean Temperature');
grid on;
hold on;
plot(year(annualT>65), annualT(annualT>65), 'ro', 'Markerfacecolor', 'r');
plot(year(annualT<60), annualT(annualT<60), 'bd', 'markerfacecolor', 'b');
warm_index = find(diff(annualT) == max(diff(annualT)));
cold_index = find(diff(annualT) == min(diff(annualT))); 
plot(year(warm_index:warm_index+1), annualT(warm_index:warm_index+1), '-r', 'linewidth', 4);
plot(year(cold_index:cold_index+1), annualT(cold_index:cold_index+1), '-b', 'linewidth', 4);
legend('Mean Temperature', 'Temperature is above 65', 'Temperature is below 60', ...
    'Fastest Temp Increase', 'Fastest Temp Decrease', 'Location', 'Northwest');
p2i = 'In general, the temperature in San Diego gets warmer over the years';
p2h = 'See figure 2';

%%Problem 3
p3a = evalc('help lottery');
p3b = lottery([2, 3, 4, 5, 6, 7]);
p3c = lottery([12, 23, 24, 34, 50, 61]);
p3d = lottery([22, 33, 44, 50, 51, 61]);
p3e = lottery([32, 43, 54, 44, 51, 61]);
p3f = lottery([42, 53, 34, 44, 51, 61]);
p3g = lottery([42, 23, 34, 44, 51, 61]);
p3h = lottery([12, 23, 34, 44, 51, 61]);

%%Problem 4
p4a = evalc('help piecewise2d');
p4b = piecewise2d(1, 1);
p4c = piecewise2d(1, -1);
p4d = piecewise2d(-1, 1);
p4e = piecewise2d(-1, -1);
p4f = piecewise2d(0, 0);
p4g = piecewise2d(0, 1);
p4h = piecewise2d(0, -1);
p4i = piecewise2d(1, 0);
p4j = piecewise2d(-1, 0);

%%Problem 5
p5a = evalc('help rgb_color');
p5b=rgb_color([1 1 1]);
p5c=rgb_color([1 0 0]);
p5d=rgb_color([0 1 0]);
p5e=rgb_color([0 0 1]);
p5f=rgb_color([1 1 0]);
p5g=rgb_color([0 1 1]);
p5h=rgb_color([1 0 1]);
p5i=rgb_color([0 0 0]);