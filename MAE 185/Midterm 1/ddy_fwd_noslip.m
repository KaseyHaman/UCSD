function dfdy = ddy_fwd_noslip(f,dy)
    
    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % forward difference
    for i=1:nx
        for j=2:ny-1
            dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
        end
    end
    
    % no slip for first and last
    dfdy(:,ny) = 0;
    dfdy(:,1) = 0;
end

