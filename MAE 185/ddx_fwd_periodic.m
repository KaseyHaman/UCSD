function [dfdx] = ddx_fwd_periodic(f, dx)
%First order Forward Difference
[a b] = size(f); %Grab size
dfdx = zeros(a, b); %preallocate
%Use for loop
for i = 1:a-1
    for j=1:b
dfdx(i, j) = (f(i+1, j) - f(i, j)) / (dx);
    end
end
%loop
for i = a
    for j = 1:b
dfdx(i, j) = (f(1, j) - f(i, j)) / (dx);
    end
end