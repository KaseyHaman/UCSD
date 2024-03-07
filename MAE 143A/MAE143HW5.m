clc;
clear all;
close all;

%%Problem 2

%Line spectra and phase
for n = 1:5
b0 = -100./(n.*pi).*(cos(2.*pi.*n)-cos(4.*pi.*n./3)) + (-150.*sin(4.*pi.*n./3)+150.*sin(2.*pi.*n./3)+100.*pi.*n.*cos(2.*pi.*n./3))/(pi^2.*n^2) ;
b(n) = abs(b0);
phase(n) = atand(b0/0);
end
n = 1:5;

%Line spectra plot
% figure(1)
% bar(n, b(n), 'BarWidth', 0.2)
% xlabel('rad/s');
% ylabel('Cn'); 
% title('Amplitude Spectrum');
% set(gca, 'FontSize', 16);
% box on; grid on;



%Phase plot
% figure(2)
% bar(n, phase(n), 'BarWidth', 0.2)
% xlabel('rad/s');
% ylabel('Phi'); 
% title('Phase');
% set(gca, 'FontSize', 16);
% box on; grid on;


%Fourier transform plot
% func = 0;
% t = -6.1:0.001:6.1;
% count = 0;
% waveform = zeros(size(t));
% a0 = 50;
% for t = -6:0.001:6
% for n=1:100
%     b0 = ( -100/(n*pi)*( 1-cos(4*pi*n/3)-cos(2*pi*n/3) ) - 150/(n.^2*pi^2)*(sin(4*pi*n/3)-sin(2*pi*n/3)) )*sin(n*pi*t/3);
%     func = func + b0;
% end
% count = count + 1;
% waveform(count) = func+ a0;
% func = 0;
% end
% 
% t = -6.1:0.001:6.1;
% figure(3)
% plot(t, waveform)
% xlabel('Time');
% ylabel('Function Value');
% title('Fourier Transform of Piecewise Function');
% xlim([-6 6])
% ylim([0 120])





% %%Problem 3
% % Fourier transform plot
% func = 0;
% t = 0:0.001:10;
% count = 0;
% waveform = zeros(size(t));
% a0 = 1/2;
% for t = 0:0.001:10
% count = count + 1;
% waveform(count) = 1/2*cos(8*pi*t)+a0;
% end
% 
% t = 0:0.001:10;
% figure(3)
% plot(t, waveform)
% xlabel('Time');
% ylabel('Function Value');
% title('Fourier Transform of Piecewise Function');
% xlim([0 2])
% ylim([-0.2 1.2])
% 
