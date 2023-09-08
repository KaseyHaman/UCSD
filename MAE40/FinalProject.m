clear all; clc;
%%Phase A
syms Vo Ic Il Ir Vs R L C Vd Vbo Ibl s
% x={Ic, Il, Ir, Vo, Vs}  <-- unknown vector
A  =[0 0 -R 1 0;      % Vo -Ir*R = 0
     1 0 0 -C*s 0;    % Ic - C*s*Vo - Vbo = 0
     0 -L*s 0 0 1/s;  % Vs/s-L*s*Il + L*Ibl = 0
     1 -1 0 0 0;      % Il = Ic 
     1 0 -1 0 0];     % Ic = Ir  
b  =[ 0; Vbo; -L*Ibl; 0; 0];
x=A\b; 

%%Phase B
syms Vo Ic Il Ir Vs R L C Vd Vbo Ibl s Vd
% x={Ic, Il, Ir, Vo, Vs}  <-- unknown vector
A  =[0 0 -R 1 0;      % Vo -Ir*R = 0
     1 0 0 -C*s 0;    % Ic - C*s*Vo - Vbo = 0
     0 -L*s 0 0 1/s;  % Vs/s-L*s*Il + L*Ibl = 0
     1 -1 -1 0 0;      % Ic = Ir + Il        
     0 0 0 1 -1/s];     %Vd/s = Vout - Vs/s
b  =[ 0; Vbo; -L*Ibl; 0; Vd/s];
x=A\b; 


%BODE Plot
% omega1 = 10;
% zeta = 5;
% P1bode = RR_tf([1, omega1, omega1.^2], [1, 2*omega1*zeta, omega1.^2]);
% figure(1), RR_bode(P1bode)



