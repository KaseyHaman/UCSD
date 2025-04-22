function [dfdx] = ddx_central(f, dx)
%Second Order Central Difference
[a b] = size(f); %Grab size
dfdx = zeros(a, b); %preallocate
%Set for loop
for i = 2:a-1
    for j = 1:b
dfdx(i, j) = (f(i+1, j) - f(i-1, j)) / (2*dx);
    end
end
%For the boundary conditions we use second order forward and backward difference
for i=1
    for j=1:b
dfdx(i, j) = (-3*f(i, j)+4*f(i+1, j)-f(i+2, j)) / (2*dx);
    end
end

for i = a
    for j=1:b
dfdx(i, j) = (3*f(i, j)-4*f(i-1, j)+f(i-2, j)) / (2*dx);
    end
end

end