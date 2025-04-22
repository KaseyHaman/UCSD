function dfdy = ddy_central_periodic(f,dy)
    
    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % central difference
    for i=1:nx
        for j=2:ny-1
            dfdy(i,j) = (f(i,j+1)-f(i,j-1))/2/dy;
        end
    end
    
    % loop
    j = 1;
    for i=1:nx
        dfdy(i,j) = (f(i,j+1)-f(i,ny))/2/dy;
    end
    
%
    j = ny;
    for i=1:nx
        dfdy(i,j) = (f(i,1)-f(i,j-1))/2/dy;
    end
    
end