function dfdy = ddy_bwd_periodic(f,dy)

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % backward difference
    for i=1:nx
        for j=2:ny
            dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
        end
    end

    % loop
    j = 1;
    for i=1:nx
        dfdy(i,j) = (f(i,j)-f(i,ny))/dy;
    end
    
end