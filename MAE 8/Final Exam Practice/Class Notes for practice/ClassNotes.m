clear all; close all; clc
%%Problem 1 HW4 Related
load('temperature.dat')
%Make figure 1 to show the average winter temperature throughout the years
year = temperature(:, 1);
winter_temp = temperature(:, [13, 2, 3]);
winter_avg_temp = mean(winter_temp');
figure(1);
plot(year, winter_avg_temp);
xlabel('year'); ylabel('winter temp (C)');
%What year is the average winter temperature coldest?
coldest_temp = min(winter_avg_temp);
p1a = year(winter_avg_temp == coldest_temp);




%%Problem 2 HW6 Related
%Use for loop to compute double/triple series. 
%Give an infinite series to estiamte pi. One of the hw problems.
%Use a nested for loop to search a string
%local min / local max / nested for-if (snow coverage)




%%Problem 3
load('survey.mat')
%Split survey into question structure
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

%Are the first and last students in the survey expecting the same letter
%grade?
p4a = isequal(survey{1}(40), survey{end}(40));
%How many students thought the midterm was moderate?
counter = 0; %my attempt
for n = 1:length(survey)
   if any(strfind(survey{n}, 'Moderate'))
       counter = counter + 1;
   else
   continue;
   end
end
counter = 0; %class answer. Correct answer is 106
for n = length(survey)
   local_ind = strfind(survey{n}, 'Moderate');
   if local_ind == 53
       counter = counter + 1;
   end
end
%Create a bar chart displaying how many students found the midterm easy,
%moderate or difficult
for i = 1:length(survey)
    survey_processed{i} = split(survey(i));
end

a = zeros(3, 1);
for i = 1:length(survey_processed)
    if strcmp(survey_processed{i}(end), 'Easy')
        a(1) = a(1)+1;
    elseif strcmp(survey_processed{i}(end), 'Moderate')
        a(2) = a(2) + 1;
    elseif strcmp(survey_processed{i}(end), 'Difficult')
        a(3) = a(3)+1;
    end
end
figure(2)
x = categorical({'Easy', 'Moderate', 'Difficult'});
bar(x, a)




%%Problem 4
load('SDweather.mat')
%Find warmest temp in 1950
for n = 1:length(SDweather)
   if SDweather(n).year == 1950
       ind_1950 = n;
       p4a = max(SDweather(n).temperature);
   end
end
%What month was the warmest temp in 1950?
p4b = find(SDweather(ind_1950).temperature == max(SDweather(ind_1950).temperature));

%What was the warmest average yearly temperature
for n = 1:length(SDweather)
    meantemp(n) = mean(SDweather(n).temperature);
end
a = max(meantemp);
b = find(meantemp == max(meantemp));
warmyear = SDweather(b).year;




%%Problem 5 Hw 7/8
% Euler or Euler Cromer method
% file input/output: importdata, fgetk, fprintf (export data)
% Compute arc length, derivative and integral





%Extra equations

%Integral = (z(2)-z(1))*(0.5*(g(1)+g(end)) + sum(g(2:end-1)));
%Arc length = sum(sqrt(diff(x).^2+diff(y).^2+diff(z).^2));
%Derivative = diff(f)./diff(x); 
%       where x = 0:0.1:10; sp2a = x; f = tanh(0.5*x).^4.*exp(-sin(x).^2);
