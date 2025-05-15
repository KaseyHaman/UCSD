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

    % forward difference for first point
    dfdy(:,1) = (f(:,1)-f(:,ny))/dy;
    
end