function [T, Z, W] = rocket(Tf, dt)
%%The function Rocket takes the duration of the flight and time step and
%%outputs the time, altitude and velocity of the rocket.
m = 10;
W(1) = 0;
Z(1) = 0;
T(1) = 0;
nstep = fix(Tf/dt);
for n = 1:nstep
  T(n+1) = T(n) + dt;
  W(n+1) = W(n) + (-gravity(Z(n)) + (thrust(T(n))/m))*dt;
  Z(n+1) = Z(n) + W(n+1)*dt;
end
end
function [Th] = thrust(T)
%%The subfunction thrust takes the time of flight of the rocket at a certain time and outputs
%%the thrust force of the rocket at that instance.
Th = 0;
if 0<=T && T<2
  Th = 670;
elseif 2<= T && T < 4
  Th = 1360;
elseif T >= 4
  Th = 0;
else
  disp('Invalid Input')
end
end
function [g] = gravity(Z)
%%The subfunction gravity takes the altitude of the rocket at a 
% certain time and outputs the force of gravity at that instance.
g = 0;
if 0 <= Z
  g = 9.81*((1-(Z/10000)^3));
elseif Z >= 10000
  g = 0;
else
  disp('Invalid Input');
end
end
