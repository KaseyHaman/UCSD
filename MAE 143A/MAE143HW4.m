clc;
clear all;
close all;

%%Problem 2
% Part e

% H = tf([2 1], [1 0])
% 
% Q = tf([2 1], [3 1])
% 
% figure(1);
% step(H)
% 
% figure(2)
% step(Q)


syms s
Y = -1/((s-3)*(s+1))
ilaplace(Y)