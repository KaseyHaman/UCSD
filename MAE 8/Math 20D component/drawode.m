function drawode(f,t_interval,t0,y0)
% drawode(f, [tmin,tmax], t0, y0)
%   
%   Plot a solution for the first order ODE y' = f(t,y)
%   using t-values from start time tmin to end time tmax
%       with initial value y(t0) = y0
%
%   f is an @ function
%     or the name of an m-file with quotes.
%   t0 need not actually be between tmin and tmax.
%
% Example: y' = (1 - t^2)*exp(y) for -4 < t < 4; y(1)=-1
%
%   f = @(t,y) (1 - t^2)*exp(y) 
%   drawode(f, [-4,4], 1, -1)

tmin=t_interval(1);
tmax=t_interval(2);
% axis manual
hold on
ax = gca;    
[numColor,~] = size(ax.ColorOrder);
if t0 ~= tmin
    [ta,ya] = ode45(f,[t0,tmin],y0);
    plot(ta,ya);
    if ax.ColorOrderIndex > 1
        ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    else
        ax.ColorOrderIndex = numColor;
    end
end
if t0 ~= tmax   
    [ta,ya] = ode45(f,[t0,tmax],y0);
    plot(ta,ya);
end
