clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';
%%Problem 1
%Preallocate
x = pi;
h1 = 1E-11;
h2 = 1E-13;
%Dh for h1 and h2
Dh1 = ((x + h1).^2 - ((x).^2))/h1;
Dh2 = ((x + h2).^2 - ((x).^2))/h2;
%Error for dh1 and dh2
Eh1 = abs(x*2-Dh1);
Eh2 = abs(x*2-Dh2);
%Q) Provide reasoning for the behavior of the errors.
 %ANS) We note that although our step size got larger, our error got
 %worse. This is likely because when we add an extremely small value to
 %x and square it, we get a lot of junk digits due to catastrophic
 %subtraction.


%%Problem 2
%Preallocate x, yprime and h
x = 1;
for n = 0:10
 h(n+1) = 10.^-n;
end
yprime = (-2*exp(-2*x)*cos(pi*x)-(exp(-2*x)*sin(pi*x)*pi));
%Solve dh, dhat and dsquiggle. Then, find the errors (e, ehat, esquiggle).
for n = 0:10
 Dh(n+1) = (exp(-2*(x+h(n+1)))*cos(pi*(x+h(n+1)))-exp(-2*x)*cos(pi*x))/h(n+1);
 Dhat(n+1) = (exp(-2*(x+h(n+1)))*cos(pi*(x+h(n+1)))-exp(-2*(x-h(n+1)))*cos(pi*(x-h(n+1))))/(2*h(n+1));
 Dsq(n+1) = ( 8* (exp(-2*(x+h(n+1)))*cos(pi*(x+h(n+1)))-exp(-2*(x-h(n+1)))*cos(pi*(x-h(n+1)))) - ...
             (exp(-2*(x+2*h(n+1)))*cos(pi*(x+2*h(n+1)))-exp(-2*(x-2*h(n+1)))*cos(pi*(x-2*h(n+1)))))/(12*h(n+1));
 eh(n+1) = abs(yprime - Dh(n+1));
 ehat(n+1) = abs(yprime - Dhat(n+1));
 esq(n+1) = abs(yprime - Dsq(n+1));
end
%Graph the errors
figure(1), hold on;
cs = 'krbgmckrbgm';
 plot(log10(h),log10(eh), 'k','LineWidth',1);
 plot(log10(h),log10(ehat), 'r','LineWidth',1);
 plot(log10(h),log10(esq), 'b','LineWidth',1);
 legend('Eh', 'Ehat', 'Esquiggle', 'Location', 'northwest')
xlabel('Log10(h)'); ylabel('log10(Error h)');
title('Log10(h) vs Log10(Error h)');
box on; grid on;
set(gca,'FontSize',10)
%Q) Where does catastrophic subtraction and taylor poly errors dominate?
  %ANS) Taylor poly is the good digitis, catastrophic subtraction errors
  %are the junk. When looking at our graph, catastrophic subtraction takes
  %over after around x = -7 for eh, x = -5 for ehat and x = -3 for
  %Esquiggle.
%Q) What does the slope tell you?
 % ANS) The larger the slope, the less error there is between that method
 % and the real function. (It also turns the functions order into slope).


%%Problem 3
%Preallocate
x = 0;
for n = 0:8
 h(n+1) = 10.^-n;
end
%Dh for h's
for n = 1:9
P3Dh(n) = (((x+h(n)).^(7/2) + x + h(n)) - (x.^(7/2) + x))/h(n);
end
%Error for dh's
for n = 1:9
P3Eh(n) = abs( ( (7/2).*(x.^(5/2)) +1 ) - P3Dh(n) );
end
%Q)Sugguest an explanation for slope behavior.
 %ans) Because our step size isn't too small, we dont get catastrophic
 %subtraction error. Because of this, as we decrease our step size, our
 %error should logically decrease. This gives us a consistent slope line.


%%Problem 4
%Preallocate x and h
x = 0;
clear h
for n = 0:5
 h(n+1) = 10.^-n;
end
ydoubleprime = -(2*(x^2-2))/((x^2+2).^2);
%Solve dh, dhat and dsquiggle. Then, find the errors (e, ehat, esquiggle).
for n = 1:6
 Dh2nd(n) = (log((x+h(n)).^2+2) + log((x-h(n)).^2+2) - 2*log((x).^2+2))/(h(n).^2);
 eh2nd(n) = abs(ydoubleprime- Dh2nd(n));
end
%Graph the errors
figure(2), hold on;
cs = 'krbgmckrbgm';
 plot(log10(h),log10(eh2nd), 'k','LineWidth',1);
 legend('eh^2', 'Location', 'northwest')
xlabel('Log10(h)'); ylabel('log10(eh^2)');
title('Log10(h) vs Log10(eh^2)');
box on; grid on;
set(gca,'FontSize',10)
%Q) Comment on result. What might the order be?
 %ans) It mirrors our other graph pretty reasonably. As for the order, its
 %a second order function as the derivative is second order. This can also
 %be proven as the slope is aprox -2.


