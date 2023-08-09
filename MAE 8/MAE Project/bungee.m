function [T, X, Y, Z, U, V, W, safety] = bungee(m, k, l, Xo, Yo, Zo, Uo, Vo, Wo)
%The function bungee takes the inputs m, k, l, Xo, Yo, Zo, Uo, Vo, Wo and
%outputs T, X, Y, Z, U, V, W, safety.
%This function utilizes a variety of physics equations to solve for unknown
%variables of the bungee jump.
global x_terrain y_terrain h_terrain Vmag
dt = 0.02;
Cd = 0.1;
A = pi;
Pa = 1.2;
g = 9.81;
h = interp2(x_terrain,y_terrain,h_terrain,Xo,Yo);
T = zeros(1, 120/0.02); T(1) = 0;
X = zeros(1, 120/0.02); X(1) = Xo;
Y = zeros(1, 120/0.02); Y(1) = Yo;
Z = zeros(1, 120/0.02); Z(1) = Zo;
U = zeros(1, 120/0.02); U(1) = Uo;
V = zeros(1, 120/0.02); V(1) = Vo;
W = zeros(1, 120/0.02); W(1) = Wo;
r = zeros(1, 120/0.02);
n = 1;
r(n) = sqrt(X(n)^2 + Y(n)^2 + Z(n)^2);
Vmag = zeros(1, 120/0.02);
Vmag(n) = sqrt(U(n)^2 + V(n)^2 + W(n)^2);
n = 0;
% Acc(n) = diff(Vmag(n));
% KE(n) = 0.5*(Vmag(n))^2;
while Z(n+1) > h
  n = n + 1;
  U(n+1) = U(n) - (k/m * (r(n)-l)/r(n) * X(n) + (Cd*Pa*A*Vmag(n))/(2*m) * U(n))*dt;
  V(n+1) = V(n) - (k/m * (r(n)-l)/r(n) * Y(n) + (Cd*Pa*A*Vmag(n))/(2*m) * V(n))*dt;
  W(n+1) = W(n) - (k/m * (r(n)-l)/r(n) * Z(n) + (Cd*Pa*A*Vmag(n))/(2*m) * W(n) + g)*dt;
  X(n+1) = X(n) + U(n+1)*dt;
  Y(n+1) = Y(n) + V(n+1)*dt;
  Z(n+1) = Z(n) + W(n+1)*dt; 
  T(n+1) = T(n) + dt;
  r(n+1) = sqrt(X(n+1)^2 + Y(n+1)^2 + Z(n+1)^2);
  Vmag(n+1) = sqrt(U(n+1)^2 + V(n+1)^2 + W(n+1)^2);
  h = interp2(x_terrain,y_terrain,h_terrain,X(n+1),Y(n+1));
if T(n+1) >= 120
  safety = 1;
  break;
else
  safety = 0;
end
end
X(X == 0) = []; Y(Y == 0) = []; Z(Z == 0) = []; Vmag(Vmag == 0) = [];
  T(end) = [];
  U(end) = [];
  V(end) = [];
  W(end) = [];
  X(end) = [];
  Y(end) = [];
  Z(end) = [];
  r(end) = [];
  Vmag(end) = [];
end

