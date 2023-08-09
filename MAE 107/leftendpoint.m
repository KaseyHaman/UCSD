function r = leftendpoint(f,a,b,n)
    h=(b-a)/n;
    % initialize steps and r.
    r=f(a);
    % need only consider the n-1 remaining sub-intervals
    for k=1:n-1
      c=a+k*h;
      r=r+f(c);
    end 
    r=h*r;
  end