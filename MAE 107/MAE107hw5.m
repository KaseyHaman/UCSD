clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';


%%Problem 3 
%Initialize functions
func1 = @(x) 0.5 + sin(pi*x) + 3*x;
func2 = @(x) 5-2*x*exp(1/(1+x^2));
%Initialize inputs
a = 0; b = 4; a1 = 0; b1 = 0; a2 = 0; b2 = 0;
e = 10^(-5);
x = (a+b)/2;
count = 1;

while abs(min(func2(x), func1(x))) >e
    if func2(x)*func2(a) > 0
        a2 = x;
        b2 = b;
    else
        b2 = x;
        a2 = a;
    end

    x2=(a2+b2)/2;

 if func1(x)*func1(a) > 0
        a1 = x;
        b1 = b;
    else
        b1 = x;
        a1 = a;
 end
    x1=(a1+b1)/2;

    if abs(func1(x1)<func2(x2))
x =(a1+b1)/2;
a = a1;
b = b1;
    else
x =(a2+b2)/2;
a = a2;
b = b2;        
    end

    count = count + 1;
    if count > 1000000
        break;
    end
end


P3Steps = count



% %%Problem 3 First Prototype
% %Initialize functions
% func1 = @(x) 0.5 + sin(pi*x) + 3*x;
% func2 = @(x) 5-2*x*exp(1/(1+x^2));
% %Initialize inputs
% a = 0; b = 4; 
% e = 10^(-5);
% x = (a+b)/2;
% count = 1;
% 
% while abs(func2(x))>e
%     if func2(x)*func2(a) > 0
%         a=x;
%     else
%         b=x;
%     end
%     x=(a+b)/2;
%     count = count + 1;
%     if count > 1000000
%         break;
%     end
%     func2steps = count;
% end





%%Problem 4 (NO CODE NEEDED, SOLVE BY HAND. CODE BELOW IS TO CHECK ANSWER.)

% func = @(x) x*exp(x) - (x^2+1);
% funcprime = @(x) x*exp(x)+exp(x)-(2*x);
% x(1) = 0;
% for k = 1:4
%     x(k+1)= x(k) - func(x(k))/funcprime(x(k));
% end



%%Problem 5
%Newtons method that converges to right solution
func = @(x) atan(x+(2*x.^3)/3) - 0.7;
funcprime = @(x) x*exp(x)+exp(x)-(2*x);
x(1) = 0;
k=1;
e = 10^(-12);
while abs(func(x(k))) > e
    x(k+1)= x(k) - func(x(k))/funcprime(x(k));
    k = k+1;
end
steps = k-1;
p5soln = (x(k))

%Viewing error
newtonerror = 0;
step = 1;
for k = 1:30
    newtonerror(k) = func(x(step+1))-func(x(step));
    step = step + 1;
end
p5error = newtonerror

%Question: Discuss errors in newtons method. 
%Answer: As we can see by our variable newtonerror, the error in our
%estimate of the root solution decreases by approxiamtely a factor of 2
%each time. This contributes to the theory of the digit doubling concept.


%Newtons method that does not converge to the right solution
clear x;
func = @(x) atan(x+(2*x.^3)/3) - 0.7;
funcprime = @(x) x*exp(x)+exp(x)-(2*x);
x(1) = 8;
k=1;
e = 10^(-12);
while abs(func(x(k))) > e
    x(k+1)= x(k) - func(x(k))/funcprime(x(k));
    k = k+1;
    if k == 1000
        break;
    end
end

p5divergeanswer = x(k)



%%Problem 6
clear x;
func = @(x) sin(x)+(5/4)*x-2;
% funcprime = @(x) cos(x)+(5/4)-2;
x(1) = 3; x(2) = 2;
funceval(1) = func(x(1)); funceval(2) = func(x(2));
%Need 7 steps to get to x8 (assuming x(1) is x0)
for n = 1:7
    x(n+2)= x(n+1) - (func(x(n+1))*(x(n+1)-x(n))) / (func(x(n+1))- func(x(n)));
    funceval(n+2) = func(x(n+2));
end

p6xvalues = x

xplot = x(1:4);
xdomain = -5:0.5:5;

%Plot the values
figure(1), hold on;
% cs = 'krbgmckrbgm';
 plot(xdomain,func(xdomain), 'k','LineWidth',1);
 plot(xplot,func(xplot), 'bo','LineWidth',1);
%  legend('Function Evalutation at x', 'Location', 'best')
xlabel('x'); ylabel('function at x');
title('Evaluating the function at x');
box on; grid on;
set(gca,'FontSize',10)

for n = 1:7
Erroranalysis(n) = funceval(n+1)/funceval(n);
end

%Q) Discuss rate of convergence
%A) When looking at the variable erroranalysis, we see that
%the error reduces in size at an increasing rate as we increase steps.
%This is due to the fact that newton's method has error of e(k+1) <=
% c*e(k)^2. Thus, our error analysis seems reasonable for the secant
% method, which is a derivation from the newton method.

%%Problem 8
%Initialize functions
e = 10^-5;
func1 = @(x) 1/2*(tan(4*x)-1/2);
func2 = @(x) atan(2*x+1/2)/4;
actualfunction1 = @(x) tan(4*x);
actualfunction2 = @(x) (2*x+1/2);
x0(1) = 0; x0(2) = 2; x0(3) = 1000;
%Loop for all x0 values
for n = 1:3

%Func 1
clear x;
x(1) = x0(n); 
x(2) = func1(x(1));
k=2;
while abs(x(k)-x(k-1)) > e
    x(k+1)= func1(x(k));
    k = k+1;
    if k == 1000
        break;
    end
end
p8func1steps(n) = k

%Func 2
clear x;
x(1) = x0(n); 
x(2) = func2(x(1));
k=2;
while abs(x(k)-x(k-1)) > e
    x(k+1)= func2(x(k));
    k = k+1;
    if k == 1000
        break;
    end
end
p8func2steps(n) = k

end


