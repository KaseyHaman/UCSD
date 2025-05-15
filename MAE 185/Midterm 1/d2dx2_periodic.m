function d2fdx2 = d2dx2_periodic(f,dx)

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    d2fdx2      = zeros(nx,ny);
    
    
    
    % central difference
    for i=2:nx-1
        for j=1:ny
            d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(i-1,j))/dx^2;
        end
    end
    
    % periodic for 1st point
    i = 1;
    for j=1:ny
        d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(nx,j))/dx^2;
    end
    
    % periodic for last point
    i = nx;
    for j=1:ny
        d2fdx2(i,j) = (f(1,j)-2*f(nx,j)+f(nx-1,j))/dx^2;
    end
    
end