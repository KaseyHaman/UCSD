function [T, X, U] = car(Tf, dt)
%%The function car translates the travelling time and time step into the time,
%%distance and velocity of the car respectively.
T(1) = 0;
X(1) = 0;
U(1) = 0;
nstep = fix(Tf/dt);
for n = 1:nstep
  T(n+1) = T(n) + dt;
  U(n+1) = U(n) + 5*sech(T(n)/20)^2*tanh(T(n)/20)*dt;
  X(n+1) = X(n) + U(n)*dt;
end
end
