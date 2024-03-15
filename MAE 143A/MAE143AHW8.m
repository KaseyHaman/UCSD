clc;
clear all;
close all;

%%Problem 3
% x = [3; 4; 1; -2]
% DFTX = fft(x)
% IDFTX = ifft(DFTX)

%%Problem 4
% A = [1; 3; 5; 0]
% B = [7; 9; 0; 0]
% conv(A, B)
% 
% 
% DFTA = fft(A)
% DFTB = fft(B)
% 
% step2 = DFTA .* DFTB
% 
% soln = ifft(step2)

%%Problem 5
% B = [1;3;5;0;0;0;0;0;0;0];
% C = [3;10;22;18;28;29;52;45;0;0];
% DFTC = fft(C)
% DFTB = fft(B)

% B = [1;4;2;6;5;3;0;0;0;0;0;0;0;0;0;0];
% C = [2;9;11;31;48;67;76;78;69;38;12;0;0;0;0;0];
% DFTC = fft(C);
% DFTB = fft(B);
% p5b = ifft(DFTC./DFTB)