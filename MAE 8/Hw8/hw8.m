clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';
hw_num = 8;

%%Problem 1

load('note.mat')
p1a = size(note);
p1b = note(:, 1);
p1c = note(end, :);
p1d = note{3, 2};
p1e = p1b;
for n = 1:length(note) %Capitlizing problem
   p1e{n}(1) = p1e{n}(1) - 32;
   a = strfind(p1e{n}, ' ');
   p1e{n}(a+1) = p1e{n}(a+1) - 32;
end

%%Problem 2

student(1).name = 'Noah Williams';
student(1).PID = 'A01';
student(1).homework = [70 91 82 93 84 85 96 78];
student(1).project = 96;
student(1).midterm = 93;
student(1).final_exam = 63;
student(2).name = 'Benjamin Frank';
student(2).PID = 'A02';
student(2).homework = [90 81 92 83 67 85 86 92];
student(2).project = 82;
student(2).midterm = 83;
student(2).final_exam = 91;
student(3).name = 'Oliver Harper';
student(3).PID = 'A03';
student(3).homework = [80 71 92 73 64 75 96 77];
student(3).project = 77;
student(3).midterm = 91;
student(3).final_exam = 76;
p2a = student(1);
p2b = student(2);
p2c = student(3);
for n = 1:3
   a = student(n).homework;
   b = find(a == min(a));
   a(b) = [];
   student(n).hw_average = mean(a);
end
p2d = student.hw_average;

%%Problem 3 
%Reading data problem

fid = fopen('class_survey.dat', 'r');
counter = 0;
   while ~feof(fid) %feof logical 0 when there's more to read
       counter = counter + 1;
       aline = fgetl(fid);
       [Q1, rest] = strtok(aline);
       [Q2, rest] = strtok(rest);
       [Q3, rest] = strtok(rest);
       [Q4, rest] = strtok(rest);
       rest = fliplr(rest);
       [Q6, Q5] = strtok(rest);
       Q5 = strtrim(fliplr(Q5));
       Q6 = strtrim(fliplr(Q6));
       survey(counter).Q1 = Q1;
       survey(counter).Q2 = Q2;
       survey(counter).Q3 = Q3;
       survey(counter).Q4 = Q4;
       survey(counter).Q5 = Q5;
       survey(counter).Q6 = Q6;
   end
   fclose(fid);
expected_grade = [survey.Q4];
%numel(strfind(expected_grade, 'A'));
target = {'A' 'B' 'C' 'D' 'f'};
for n = 1:length(target)
   number_student(n) = numel(strfind(expected_grade,target{n}));
end
figure(1)
bar(number_student);
xlabel('Letter Grade')
ylabel('Number of Students')
set(gca, 'XTickLabel', target);
title('Expected Course Grade')

coding_experience = [survey.Q2];
p3b = numel(strfind(coding_experience, 'Yes'));

study_hours = [survey.Q5];
a = numel(strfind(study_hours, '7 - 9'));
b = numel(strfind(study_hours, '> = 10'));
p3c = a + b;

% class_level = [survey.Q1];
% sophomore_index = strfind(class_level, 'Sophomore');


%sophomores that expect an A
p3d = 0;
sophomore_string = 'Sophomore';
for n = 1:155
    if any(strfind(survey(n).Q1, sophomore_string))
        p3d = p3d + numel(strfind(survey(n).Q4, 'A'));
    end
end

      


