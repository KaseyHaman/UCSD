function d2fdx2 = d2dx2periodic(f,dx)

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    d2fdx2      = zeros(nx,ny);
    
    dx2         = dx^2;
    
    % central difference
    for i=2:nx-1
        for j=1:ny
            d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(i-1,j))/dx2;
        end
    end
    
    % Boundary periodic
    i = 1;
    for j=1:ny
        d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(nx,j))/dx2;
    end
    
%
    i = nx;
    for j=1:ny
        d2fdx2(i,j) = (f(1,j)-2*f(i,j)+f(i-1,j))/dx2;
    end
    
end