function r = trap(f,a,b,n)
    % initialize h and r
    h=(b-a)/n;
    r=0;
    % need only consider the n-1 remaining sub-intervals
    for k=1:n-1
      c=a+k*h;
      r=r+f(c);
    end 
    %Formula for trap
    r=h*((f(a)+f(b))/2 +r);

  end