clc;
clear all;
close all;

%%Problem 1
% 
% syms s
% F = tf([1 0], [1 10^4])
% bode(F)

%%Problem 3
syms s
% F = tf([-0.01 0], [-20 -1000 -9.81*10^-2])
% F = (-s/(-20*s^2-1000*s-980))
% ilaplace(F)

syms t b
B = (t*cos(b*t))
laplace(B)