clear all;
clc
RRbase='/Users/kasey/RR'; cd(RRbase); RR_path_init
cd '/Users/kasey/Desktop/School Related/UCSD/MAE40'


%%Problem 1
syms s Vi Vs C L R R1 R2
% x={I_L, I_c, I_R, Vo}  <-- unknown vector
A  =[ -1 1 1 0;     % I_L - Ic - Ir = 0
     L*s  0  0 1;   % L*s*I_L + Vo = Vi
      0  -1 0 C*s;  % C*s*Vo - Ic = 0
      0 0 R -1];    % Ir * R - Vo = 0
b  =[ 0; Vi; 0; 0];
x=A\b; Vo_LPF2_P1=simplify(x(4))

%%Problem 2
syms zeta
zeta = 0.1;
omega4=10; F_LPF2_P2=RR_tf([omega4^2],[1 2*zeta*omega4 omega4^2]);
figure(1), RR_bode(F_LPF2_P2)

zeta = 0.7;
omega4=10; F_LPF2_P2=RR_tf([omega4^2],[1 2*zeta*omega4 omega4^2]);
figure(2), RR_bode(F_LPF2_P2)

zeta = 1;
omega4=10; F_LPF2_P2=RR_tf([omega4^2],[1 2*zeta*omega4 omega4^2]);
figure(3), RR_bode(F_LPF2_P2)

%%Problem 3
syms s V1 V2 V3 C L Rd Cd
% x={I_L, I_c, I_d, V2, V3}  <-- unknown vector
A  =[ 1 -1 1 0 0;       % I_L - Ic - Id = 0
     0 s*L 0 1 0;       % L*s*I_L + V2 = V1
      -1 0 0 C*s 0;     % C*s*V2 - Ic = 0
      0 0 Rd -1 1;      % I_d * Rd + V2 - V3 = 0
      0 0 -1 0 Cd*s];   % Cd*s*V3 - I_d = 0
b  =[ 0; V1; 0; 0; 0];
x=A\b; Vo_LPF2_P1=simplify(x(4))

%%Problem 4
syms s V1 V2 V3 C L Rd Cd
% x={I_L, I_c, I_d, V2, V3}  <-- unknown vector
A  =[ 1 -1 1 0 0;            % I_L - Ic - Id = 0
     0 s*L 0 1 0;            % L*s*I_L + V2 = V1
      -1 0 0 C*s 0;          % C*s*V2 - Ic = 0
      0 0 sqrt(1/L*C) -1 1;  % I_d * Rd + V2 - V3 = 0
      0 0 -1 0 4*C*s];       % Cd*s*V3 - I_d = 0
b  =[ 0; V1; 0; 0; 0];
x=A\b; Vo_LPF2_P1=simplify(x(4))
