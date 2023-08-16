clear all;
clc
RRbase='/Users/kasey/RR'; cd(RRbase); RR_path_init
cd '/Users/kasey/Desktop/School Related/UCSD/MAE40'

%Note, the passive filters help simplify the transfer function, RR_tf
%prepares it for bode plot.

% RR_Ex10_02_passive_filters


%%PRACTICE taking BODE plot of a transfer function
syms zeta
zeta =0.1;
%take any zeta 0<1 (idk if its necessarily 0 to 1 but thats the class ex)
G1=RR_tf([1 10 10000], [100 (200*zeta+1) (2*zeta+100) 1])
%G1 RR_tf ex.... s^2 + 10s + 10000 as numerator.
RR_bode(G1)


%%Example for solving and graphing G(s) = w1 / (s + w1)
omega1 = 10;
test = RR_tf([omega1], [1, omega1]);
RR_bode(test)

