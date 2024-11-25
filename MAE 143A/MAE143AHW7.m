timedd = [1:1:200]
for n = 1:1:200
    angled(n) = angled(n)-180
end
figure(1)
plot(timedd, angled)
figure(2)
plot(timedd, abs(THEFFT))

% clc;
% clear all;
% close all;

% %%Problem 1
% p1a = tf([1 0 0], [1 -0.2 -0.24], -1);
% figure (1)
% bode(p1a)

THEFFT = fft(angled)