clc;
clear all;
close all;


syms s
F = 3*(s-0.9)/(s*(s+1.2)*(s+0.9)*(s+1))


%%Problem 2

S = tf('s');
SYS1 = tf([4 5],[1 5 4]);
SYS2 = tf([4 0],[1 0 5]);
h1 = @(t) (11/3).*exp(-4.*t) + (1/3).*exp(-t);
h2 = @(t) 4.*cos(sqrt(5).*t);
y1 = @(t) 5/4 - (11/12).*exp(-4.*t) - (1/3).*exp(-t);
y2 = @(t) (4/sqrt(5)).*sin(sqrt(5).*t);
t = 0.000001:0.00001:3;
%
figure (1)
title('Calculated vs. Actual Responses')
subplot(1,2,1)
hold on;
impulseplot(SYS1,SYS2,4)
plot(t,h1(t),'color','red','LineWidth',1.5)
plot(t,h2(t),'color','cyan','LineWidth',1.5)
legend('Impulse Response Actual Part I','Impulse Response Actual Part II', ...
   'Impulse Response Calculated Part I','Impulse Response Calculated Part II', 'location', 'southwest')
subplot(1,2,2)
hold on;
stepplot(SYS1,SYS2,4)
plot(t,y1(t),'color','red','LineWidth',1.5)
plot(t,y2(t),'color','cyan','LineWidth',1.5)
legend('Actual Step Response Part I','Actual Step Response Part II', ...
   'Calculated Step Response Part I','Calculated Step Response Part II', 'location', 'southwest')



%Problem 4

%G(jw) plot
w = 0:0.001:10;
G_mag = 1./(sqrt(w.^2 + 1));
G_phase = -atand(w);

figure(2)
subplot(2,1,1)
plot(w, G_mag, 'linewidth', 2);
grid;
xlabel('\omega (rad/s)');
ylabel('|G(j\omega)|');
title('Magnitude of G(s)');

subplot(2,1,2)
plot(w, G_phase, 'linewidth', 2);
grid;
xlabel('\omega (rad/s)');
ylabel('\angle G(j\omega)');
title('Phase of G(s)');

%G(e^jwTs) plot

Ts = 0.01; %0.01 seconds

w = 0:pi/1000:pi;

b = [Ts, Ts];
a = [2+Ts, Ts-2];

G = freqz(b,a,w);
G_mag = abs(G);
G_phase = angle(G)*180/pi;

figure(3)
subplot(2,1,1)
plot(w/pi,G_mag);
grid
xlabel('Frequency \omega (\times \pi (rad/sample)');
ylabel('Magnitude');
title('Magnitude Response of G_d(z)  T_s = 0.01');

subplot(2,1,2)
plot(w/pi, G_phase);
grid
xlabel('Frequency \omega (\times \pi (rad/sample)');
ylabel('Phase (degrees)');
title('Phase Response of G_d(z)   T_s = 0.01');



Ts = 0.1; %0.1 seconds

w = 0:pi/1000:pi;

b = [Ts, Ts];
a = [2+Ts, Ts-2];

G = freqz(b,a,w);
G_mag = abs(G);
G_phase = angle(G)*180/pi;

figure(4)
subplot(2,1,1)
plot(w/pi,G_mag);
grid
xlabel('Frequency \omega (\times \pi (rad/sample)');
ylabel('Magnitude');
title('Magnitude Response of G_d(z)  T_s = 0.1');

subplot(2,1,2)
plot(w/pi, G_phase);
grid
xlabel('Frequency \omega (\times \pi (rad/sample)');
ylabel('Phase (degrees)');
title('Phase Response of G_d(z)   T_s = 0.1');


Ts = 1; %1 second

w = 0:pi/1000:pi;

b = [Ts, Ts];
a = [2+Ts, Ts-2];

G = freqz(b,a,w);
G_mag = abs(G);
G_phase = angle(G)*180/pi;

figure(5)
subplot(2,1,1)
plot(w/pi,G_mag);
grid
xlabel('Frequency \omega (\times \pi (rad/sample)');
ylabel('Magnitude');
title('Magnitude Response of G_d(z)  T_s = 1');

subplot(2,1,2)
plot(w/pi, G_phase);
grid
xlabel('Frequency \omega (\times \pi (rad/sample)');
ylabel('Phase (degrees)');
title('Phase Response of G_d(z)   T_s = 1');


%Problem 5
%Problem i)
sys = tf([0 3 -3.6],[1 -1.9 0.9], Ts=-1);
t = 0:1:100; u = exp(-t); 
figure(6)
lsim(sys,u, t)  
%Problem ii)
sys = tf([0 3 2.7],[1 -2.1 1.08 0], Ts=-1);
t = 0:1:100; u = exp(-t); 
figure(7)
lsim(sys,u, t) 
%Problem iii)
sys = tf([0 3 -2.7],[1 2.1 1.08 0], Ts=-1);
t = 0:1:100; u = exp(-t); 
figure(8)
lsim(sys,u, t)
