clear all; close all; clc; format long
name = 'Kasey Haman';
id = 'A16978114';
hw_num = 'final';
form = 'B';

%% Problem 1: 
load('temperature.dat')
 year1992 = find(temperature(:, 1) == 1992);
 p1a = max(temperature(year1992, 2:13));
 p1b = find(temperature(year1992, 2:13) == p1a);
 maxtemp = max(max(temperature(:, 2:13)));
 p1c = maxtemp;
[r, c] = find(temperature(:, 2:13)==max(max(temperature(:, 2:13))));
p1d = temperature(r, 1);

year2012 = find(temperature(:, 1) == 2012);
p1e = numel(find(temperature(year1992:year2012, 2:13) > 71));

year = temperature(:, 1);
autumn_temp = temperature(:, [10:12]);
autumn_avg_temp = mean(autumn_temp');
p1f = autumn_avg_temp;

figure(1);
plot(year, autumn_avg_temp); hold on;
xlabel('Year'); ylabel('Autumn temp (C)');
autumnhigh = max(autumn_avg_temp);
autumnlow = min(autumn_avg_temp);
[highyear] = find(autumn_avg_temp == autumnhigh);
[lowyear] = find(autumn_avg_temp == autumnlow);
plot(temperature(highyear, 1), autumn_avg_temp(highyear),'ro', 'markerfacecolor', 'r', 'markersize', 10)
plot(temperature(lowyear, 1), autumn_avg_temp(lowyear),'bo', 'markerfacecolor', 'b', 'markersize', 10)
xlabel('Years'); ylabel('Temperature (f)');
title('Average Autumn Temperature');
legend('Average Temperature', 'Highest Average Temperature', 'Lowest Average Temperature', 'Location', 'best')
axis([1852 2020 58 74])

p1g = 'See figure 1';
%% Problem 2: 
p1a = 0;
for k = 1:40
    for l = 1:50
        for m = 1:60
            p1a = p1a + 1/(2^k + 2^l + 2^m);
        end
    end
end

cvpos = pi + (10^(-6)); 
cvneg = pi - (10^(-6));

p2b = zeros(1, 200);
p2b(1) = 1;
k = 200;
for n = 1:k
    p2b(n+1) = p2b(n) + (2^n*factorial(n)^2)/factorial(2*n+1);
end
p2b(:) = p2b(:) * 2;
a = find(cvneg< p2b & p2b< cvpos);
p2b = p2b(a(1));
p2c = a(1) - 1;

x = [-11:0.1:11];
y = [-12:0.1:12];
for n = 1:length(x)
     for m = 1:length(y)
        f(n, m) = exp(-(cos(x(n)./2) + sin(y(m)./3))^2);
    end
end

figure(2)
surf(x, y, f'); shading interp; hold on; view(3)
counter_lmax = 0;
for n = 2:length(x)-1
     for m = 2:length(y)-1
          nb = f(n-1:n+1, m-1:m+1);
          nb = nb(:);
         if f(n, m) >= max(nb)
             counter_lmax = counter_lmax + 1;
             x_localmax(counter_lmax) = x(n);
             y_localmax(counter_lmax) = y(m);
             f_localmax(counter_lmax) = f(n, m);
         end
     end
end
plot3(x_localmax, y_localmax, f_localmax, 'bo', 'Markersize', 5, 'MarkerFacecolor', 'b');
legend('Terrain', 'Local Maximums', 'location', 'best')
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Terrain of f');


%% Problem 3: 
load('survey.mat')
counter = 0;
for n = 1:size(survey)
    counter = counter + 1;
       [Q1, rest] = strtok(survey{n});
       [Q2, rest] = strtok(rest);
       [Q3, rest] = strtok(rest);
       [Q4, rest] = strtok(rest);
       rest = fliplr(rest);
       [Q6, Q5] = strtok(rest);
       Q5 = strtrim(fliplr(Q5));
       Q6 = strtrim(fliplr(Q6));
       s(counter).Q1 = Q1;
       s(counter).Q2 = Q2;
       s(counter).Q3 = Q3;
       s(counter).Q4 = Q4;
       s(counter).Q5 = Q5;
       s(counter).Q6 = Q6;
end

a = strfind(s(2).Q1, 'Freshman');
b = strfind(s(end).Q1, 'Freshman');
p3a = (a == b);
a = strfind(s(16).Q1, 'Freshman');
b = strfind(s(146).Q1, 'Freshman');
p3b = (a == b);

difficulty =[s.Q6];
p3c = numel(strfind(difficulty, 'Moderate'));

p3d = 0;
freshman_string = 'Freshman';
for n = 1:155
    if any(strfind(s(n).Q1, freshman_string))
        p3d = p3d + numel(strfind(s(n).Q6, 'Moderate'));
    end
end

class_level = [s.Q1];
target = {'Freshman' 'Sophomore' 'Junior' 'Senior' 'Null'};
for n = 1:length(target)
   number_student(n) = numel(strfind(class_level,target{n}));
end
figure(3)
bar(number_student);
ylabel('Number of Students')
set(gca, 'XTickLabel', target);
title('Class Diversity')

%% Problem 4: 
load('SDweather.mat')
for n = 1:length(SDweather)
   if SDweather(n).year == 1970
       ind_1970 = n;
       p4a = min(SDweather(n).temperature);
   end
end
p4b = find(SDweather(ind_1970).temperature == min(SDweather(ind_1970).temperature));


for n = 1:length(SDweather)
    meanrainfall(n) = (SDweather(n).annual_rainfall_avg);
end
a = max(meanrainfall);
p4c = a;
b = find(meanrainfall == max(meanrainfall));
p4d = SDweather(b).year;
minvalue = min(meanrainfall);
c = find(meanrainfall == min(meanrainfall));
minyear = SDweather(c).year;

rainyears = [SDweather.year];
rainamount = [SDweather.annual_rainfall_avg];
figure(4)
bar(rainyears, rainamount); hold on;
plot(p4d, p4c, 'md', 'markerfacecolor', 'm', 'markersize', 10)
plot(minyear, minvalue, 'gd', 'markerfacecolor', 'g', 'markersize', 10)
title('Annual Average Rainfall')
xlabel('Year'); ylabel('Average Rainfall (inch)');

p4e = 'See figure 4';

%% Problem 5: 

c = [0 0.6 1.2 1.8 2.4 3.2];
for n = 1:length(c)
    [T{n}, X{n}, V{n}] = spring_mass_damper(c(n));
end

figure(5); hold on;
cs = 'krbgmc';
for n = 1:length(c)
    plot(T{n},X{n},cs(n),'LineWidth',2);    
    legend_string{n} = sprintf('c = %3.1f',c(n));
end
xlabel('Time (s)'); ylabel('Displacement (m)'); 
legend(legend_string,'Location','Northeast');
title('Effect of Spring Mass Damper');
box on; grid on;
set(gca,'FontSize',14)
p5a = 'See figure 5';

x_tip{1} = 0;
t_tip{1} = 0;
for m = 1:length(c)
    counter = 0;
for n = 2:length(X{m})-1
    if X{m}(n) > (X{m}(n+1)) && X{m}(n) > X{m}(n-1)
        counter = counter  + 1;
        x_tip{m}(counter) = X{m}(n);
        t_tip{m}(counter) = T{m}(n);
    end
end
end

p5b{1} = 0;
for m = 1:length(c)
    p5b{m} = 2*pi/(t_tip{m}(2) - t_tip{m}(1));
end
    
p5b = [p5b{:}];
