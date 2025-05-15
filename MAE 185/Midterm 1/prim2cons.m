function U = prim2cons(rho, u, v,T, cv)
%Function for converting primitive variable to conservative variable
%Returns single array 'U' of size [4,nx,ny] that contains the conservative variables

[nx,ny] = size(rho);

U = zeros(4,nx,ny);

%Rho
U(1,:,:) = rho;

%Rho u
U(2,:,:) = rho.*u;

%Rho v
U(3,:,:) = rho.*v;

%Total energy
U(4,:,:) = rho.*(cv*T+(u.^2+v.^2)/2);

end

