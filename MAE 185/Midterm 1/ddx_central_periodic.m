function dfdx = ddx_central_periodic(f,dx)
% determine field size
[nx,ny] = size(f);
% allocate return field
dfdx = zeros(nx,ny);
% central difference
    for i=2:nx-1
        for j=1:ny
            dfdx(i,j) = (f(i+1,j)-f(i-1,j))/2/dx;
        end
    end
% forward difference for first point
i = 1;
    for j=1:ny
        dfdx(i,j) = (f(i+1,j)-f(nx,j))/2/dx;
    end
% backward difference for last point
i = nx;
    for j=1:ny
        dfdx(i,j) = (f(1,j)-f(i-1,j))/2/dx;
    end
end

