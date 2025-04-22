function [dfdx] = ddx_bwd_periodic(f, dx)
%First order Backward Difference
[a b] = size(f); %Grab size
dfdx = zeros(a, b); %preallocate
%Use for loop

for i=2:a
    for j=1:b
dfdx(i, j) = (f(i, j) - f(i-1, j)) / (dx);
    end
end
%For the boundary condition we loop
for i =1
    for j = 1:b
dfdx(i, j) = (f(i, j) - f(a, j)) / (dx);
    end
end

end