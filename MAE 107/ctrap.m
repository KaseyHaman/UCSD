function r = ctrap(f,a,b,n)
 % initialize h and r
    h=(b-a)/n;
    r=0;
    % need only consider the n-1 remaining sub-intervals
    for k=1:n-1
      c=a+k*h;
      r=r+f(c);
    end 
    %Formula for ctrap
    r=h*( (f(a)+f(b)) /2 + r - 1/24*( 3*f(a)-4*f(a+h)+f(a+2*h)+f(b-2*h)-4*f(b-h)+3*f(b) ));

  end