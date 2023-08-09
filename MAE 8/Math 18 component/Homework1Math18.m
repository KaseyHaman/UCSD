clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
% 
% a = 1; e = 5; h = 8; k = 11; m = 13; n = 14; s = 19; y = 25;
% kaseyhaman = k + a + s + e + y + h + a + m + a + n;

f = @(x,y) x^2 - 1;
h = 0.2;
[x,y] = Euler(h, 0, 1, 2, f);

plot(x, y, 'b')
hold on
