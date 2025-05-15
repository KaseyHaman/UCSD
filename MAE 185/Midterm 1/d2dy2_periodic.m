function d2fdy2 = d2dy2_periodic(f,dy)
[nx,ny] = size(f);
% allocate return field
d2fdy2 = zeros(nx,ny);
%central difference
    for i=1:nx
        for j=2:ny-1
            d2fdy2(i,j) = (f(i,j+1)-2*f(i,j) + f(i,j-1))/dy^2;
        end
    end
%periodic for 1st point
j = 1;
    for i=1:nx
        d2fdy2(i,j) = (f(i,j+1)-2*f(i,j)+f(i,ny))/dy^2;
    end
%periodic for last point
j = ny;
for i = 1:nx
    d2fdy2(i,j) = (f(i,1)-2*f(i,j)+f(i,j-1))/dy^2;
end
end