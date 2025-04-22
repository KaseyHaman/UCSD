function [dfdx] = ddx_central_periodic(f, dx)
%Second Order Central Difference
[a b] = size(f); %Grab size
dfdx = zeros(a, b); %preallocate
%Set for loop
for i = 2:a-1
    for j = 1:b
dfdx(i, j) = (f(i+1, j) - f(i-1, j)) / (2*dx);
    end
end
%loop
for i=1
    for j=1:b
dfdx(i, j) = (f(i+1, j) - f(a, j)) / (2*dx);
    end
end
%
for i = a
    for j=1:b
dfdx(i, j) = (f(1, j) - f(i-1, j)) / (2*dx);
    end
end

end