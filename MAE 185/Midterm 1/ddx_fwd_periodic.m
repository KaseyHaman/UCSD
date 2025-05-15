%Function for calculating first order forward difference
%Inputs: f - An array containing all grid points
%dx - grid spacing

function dfdx = ddx_fwd_periodic(f,dx)
% determine field size
[nx,ny] = size(f);
% allocate return field
dfdx = zeros(nx,ny);
% forward difference
    for i=1:nx-1
        for j=1:ny
            dfdx(i,j) = (f(i+1,j)-f(i,j))/dx;
        end
    end
% periodic BC for last point
dfdx(nx,:) = (f(nx,:)-f(1,:))/dx;