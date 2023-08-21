clear all;
clc
RRbase='/Users/kasey/RR'; cd(RRbase); RR_path_init
cd '/Users/kasey/Desktop/School Related/UCSD/MAE40'

%%Problem 1
syms s Vo V1 V2 Cc Cb Rd Ra
% x={I_o, I_a, I_b, I_c, I_d, V1 V2}  <-- unknown vector
A  =[ 1 -1 -1 0 0 0 0;           % Io = Ia + Ib 
     1 -1 0 -1 0 0 0;           % Io = Ia + Ic
     0 0 1 -1 -1 0 0;          % Ib = Id + Ic
     0 Ra 0 0 0 1 0;           % Vo = V1 + Ia *Ra
     0 0 1/(Cb*s) 0 0 0 1;     % Ib / (Cb * s) + V2 = Vo
     0 0 0 0 -Rd 0 1;          % V2 - Id * Rd = 0
     0 0 0 1 0 Cc*s -Cc*s];    % Ic - Cc*s(V2 - V1) = 0
b  =[ 0; 0; Vo; Vo; 0; 0; 0];
x=A\b; Vo_LPF2_P1=simplify(x(6))