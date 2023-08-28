clear all;
clc

%%Problem 1
%Method 1
syms s R Rload L C c1 Vo           % Laplace variable s, parameters, input V0
syms Vi Vm Ir Irload Il Ic Iload    % variables to be solved for (output is V1)
eqn1= Il == (Vi-Vm)/(L*s)        % Component eqns
eqn2= Ic == C*s*(Vm-Vo)
eqn3=    Ir == Vo/R
eqn4=    Irload == Vo/Rload
eqn5=    Il == Ic    % KCLs
eqn6=    Ic == Ir + Irload
A=solve(eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,Vo,Vm,Ir,Irload,Ic,Il); % do the solve
A.Vo/Vi    % output of interest (V0), divided by input (V1)

%Method 2
syms s Vi C Rload R L
% x=
A  =[L*s 0 0 0 1 0;        
     0 -1 0 0 C*s -C*s;         
     0 0 -1 0 0 1/R;        
     0 0 0 -1 0 1/Rload;           
     1 -1 0 0 0 0;     
     0 -1 1 1 0 0];   
b  =[ Vi; 0; 0; 0; 0; 0];
x=A\b; Vo_LPF2_P1=simplify(x(6))