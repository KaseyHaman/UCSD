clc;
clear all;
close all;

%%Problem 1
p1a = tf([1 0 0], [1 -0.2 -0.24]);
figure (1)
bode(p1a)