function f = piecewise2d(x, y)
%Automatic piecewise function
%   The function sorts x and y values into their respective function (based
%   on their domain values). Ex. piecewise2d(1, 1) will run the function 
%   f = (10*x + 10*y) while piecewise2d(1, -1) will run the function 
%   f = (10*x -10*y).
if x >= 0  && y > 0
    f = (10*x + 10*y);
 
elseif x < 0 && y >= 0
    f = (-10*x + 10*y);
 
elseif x <= 0 && y < 0
    f = (-10*x -10*y);
 
elseif x > 0 && y <= 0
    f = (10*x -10*y);
else 
    f = 0;
end 

end