function [dfdx] = ddx_fwd(f, dx)
%First order Forward Difference
[a b] = size(f); %Grab size
dfdx = zeros(a, b); %preallocate
%Use for loop
for i = 1:a-1
    for j=1:b
dfdx(i, j) = (f(i+1, j) - f(i, j)) / (dx);
    end
end
%For the boundary condition we use first order backward difference
for i = a
    for j = 1:b
dfdx(i, j) = (f(i, j) - f(i-1, j)) / (dx);
    end
end
end