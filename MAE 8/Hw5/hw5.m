clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 5;

%%Problem 1
p1a = evalc('help assign_grade');
load('studentA.mat'); [p1b, p1c] = assign_grade(homework, midterm, project, final);
load('studentB.mat'); [p1d, p1e] = assign_grade(homework, midterm, project, final);
load('studentC.mat'); [p1f, p1g] = assign_grade(homework, midterm, project, final);

%%Problem 2
p2a=evalc('help days_in_month');
p2b=days_in_month('jan',0);
p2c=days_in_month('feb',0);
p2d=days_in_month('feb',1);
p2e=days_in_month('apr',0);
p2f=days_in_month('aug',1);
p2g=days_in_month('oct',0);
p2h=days_in_month('nov',1);
p2i=days_in_month('Dec',0);

%%Problem 3
n = 49;
a = 0;
b = 1;
p3a = 1:50;
for i = 1:n
    c = a+b;
    a = b;
    b = c;
    p3a(i+1) = c;
end

%sum(p3a);
p3b = 0;
for n = 1:50
    p3b = p3b + p3a(n);
end

n = 49;
p3c = 1:50;
p3c(1) = 0;
for i = 1:n
p3c(i+1) = (p3a(i+1)./p3a(i));
end

%%Problem 4
n = 10;
p4a = 1:10;
p4a(1) = sqrt(2);
for i = 2:n
     p4a(i) = sqrt(2 + (p4a(i-1)));
end
p4a = p4a(end);

n = 9;
p4b = 1:9;
p4b(1) = sqrt(1+2);
for i = 2:n
     p4b(i) = sqrt(1+(2.*(p4b(i-1))));
end
p4b = p4b(end);

n = 10;
p4c = 1:10;
p4c(1) = sqrt(2);
for i = 2:n
    if rem(i, 2) == 0
        p4c(i) = sqrt(2 - p4c(i-1));
    else
        p4c(i) = sqrt(2 + p4c(i-1));
    end
end
p4c = p4c(end);

%%Problem 5
p5 = evalc('type("survey.dat")');
