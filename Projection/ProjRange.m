function y = ProjRange(x,m,M) ; 
    y = x ; 
    y(x>M) = M ; 
    y(x<m) = m ; 
end