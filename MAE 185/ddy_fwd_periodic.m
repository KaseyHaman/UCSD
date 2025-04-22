function dfdy = ddy_fwd_periodic(f,dy)

     % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % forward difference
    for i=1:nx
        for j=1:ny-1
            dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
        end
    end
    
    % loop
    j = ny;
    for i=1:nx
        dfdy(i,j) = (f(i,1)-f(i,j))/dy;
    end
    
end