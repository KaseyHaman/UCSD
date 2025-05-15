function U = prim2cons(rho, u, v, T, cv)
%Converts primitive variables to conservative variables

%Define e and solve for Et
e = T.*cv;
Et = rho.*(e + (u.^2 + v.^2)./2); 
 
%Sort U
U = zeros(4, size(rho, 1), size(rho, 2));
U(1, :, :) = rho;
U(2, :, :) = rho.*u;
U(3, :, :) = rho.*v;
U(4, :, :) = Et;

end