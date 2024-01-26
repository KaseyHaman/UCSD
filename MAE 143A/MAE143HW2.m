clc;
clear all;
close all;

%%Problem 2

% %Part A
% sys_a = tf([0 4 5], [1, 5, 4], Ts=-1)
% impulse(sys_a)
% step(sys_a)
% 
% %Part B
% sys_b = tf([0 4 0], [1 0 5], Ts=-1)
% impulse(sys_b)
% step(sys_b)
% 


%%Problem 4
% 
% %G(jw) plot
% w = 0:0.001:10;
% G_mag = 1./(sqrt(w.^2 + 1));
% G_phase = -atand(w);
% 
% subplot(2,1,1)
% plot(w, G_mag, 'linewidth', 2);
% grid;
% xlabel('\omega (rad/s)');
% ylabel('|G(j\omega)|');
% title('Magnitude response of G(s)');
% 
% subplot(2,1,2)
% plot(w, G_phase, 'linewidth', 2);
% grid;
% xlabel('\omega (rad/s)');
% ylabel('\angle G(j\omega)');
% title('Phase response of G(s)');

% %G(e^jwTs) plot
% 
% Ts = 0.01; %0.01 seconds
% 
% w = 0:pi/1000:pi;
% 
% b = [Ts, Ts];
% a = [2+Ts, Ts-2];
% 
% G = freqz(b,a,w);
% G_mag = abs(G);
% G_phase = angle(G)*180/pi;
% 
% subplot(2,1,1)
% plot(w/pi,G_mag);
% grid
% xlabel('Normalized Frequecny \omega (\times \pi (rad/sample)');
% ylabel('Magnitude');
% title('Magnitude Response of G_d(z)  T_s = 0.01');
% 
% subplot(2,1,2)
% plot(w/pi, G_phase);
% grid
% xlabel('Normalized Frequecny \omega (\times \pi (rad/sample)');
% ylabel('Phase (degrees)');
% title('Phase Response of G_d(z)   T_s = 0.01');
% 
% 
% 
% Ts = 0.1; %0.1 seconds
% 
% w = 0:pi/1000:pi;
% 
% b = [Ts, Ts];
% a = [2+Ts, Ts-2];
% 
% G = freqz(b,a,w);
% G_mag = abs(G);
% G_phase = angle(G)*180/pi;
% 
% subplot(2,1,1)
% plot(w/pi,G_mag);
% grid
% xlabel('Normalized Frequecny \omega (\times \pi (rad/sample)');
% ylabel('Magnitude');
% title('Magnitude Response of G_d(z)  T_s = 0.1');
% 
% subplot(2,1,2)
% plot(w/pi, G_phase);
% grid
% xlabel('Normalized Frequecny \omega (\times \pi (rad/sample)');
% ylabel('Phase (degrees)');
% title('Phase Response of G_d(z)   T_s = 0.1');
% 
% 
% Ts = 1; %1 second
% 
% w = 0:pi/1000:pi;
% 
% b = [Ts, Ts];
% a = [2+Ts, Ts-2];
% 
% G = freqz(b,a,w);
% G_mag = abs(G);
% G_phase = angle(G)*180/pi;
% 
% subplot(2,1,1)
% plot(w/pi,G_mag);
% grid
% xlabel('Normalized Frequecny \omega (\times \pi (rad/sample)');
% ylabel('Magnitude');
% title('Magnitude Response of G_d(z)  T_s = 1');
% 
% subplot(2,1,2)
% plot(w/pi, G_phase);
% grid
% xlabel('Normalized Frequecny \omega (\times \pi (rad/sample)');
% ylabel('Phase (degrees)');
% title('Phase Response of G_d(z)   T_s = 1');


%%Problem 5
% 
% sys = tf([0 3 -3.6],[1 -1.9 0.9], Ts=-1)
% t = 1;
% lsim(sys, ([1; 1]), t)

