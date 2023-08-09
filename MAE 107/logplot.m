function logplot(a, n, errorLE, errorTrap, errorCTrap)
%This function plots the errors of our estimated integrals.

%Set problem # for string title
exp_num = a + 4;
strtitle = sprintf('Problem #%1i, Error of Integral Approximation', exp_num);

if errorCTrap == 0 %Purpose of if/else is to fix the legend when Ctrap = 0 for problem 6.

figure(a), hold on;
cs = 'krbgmckrbgm';
 plot(log10(n),log10(errorLE), '-ko','LineWidth',1);
 plot(log10(n),log10(errorTrap), 'r','LineWidth',1);
 plot(log10(n),log10(errorCTrap), 'b','LineWidth',1);
 legend('Left Endpoint', 'Trapezoid', 'Location', 'southwest')
xlabel('Log10(n)'); ylabel('log10(Error)');
title(strtitle);
box on; grid on;
set(gca,'FontSize',10)

else
    
figure(a), hold on;
cs = 'krbgmckrbgm';
 plot(log10(n),log10(errorLE), '-ko','LineWidth',1);
 plot(log10(n),log10(errorTrap), 'r','LineWidth',1);
 plot(log10(n),log10(errorCTrap), 'b','LineWidth',1);
 legend('Left Endpoint', 'Trapezoid', 'Corrected Trapezoid', 'Location', 'southwest')
xlabel('Log10(n)'); ylabel('log10(Error)');
title(strtitle);
box on; grid on;
set(gca,'FontSize',10)


end


end