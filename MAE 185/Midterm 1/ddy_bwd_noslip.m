function dfdy = ddy_bwd_noslip(f,dy)

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % backward difference
    for i=1:nx
        for j=2:ny-1
            dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
        end
    end

    %no slip for first and last
    dfdy(:,1) = 0;
    dfdy(:,ny) = 0;
    
end