function C= cellgen( f, n )
%cellgen(f,n) Create a cell C of size n, where C{i}=f(i); 

C=cell(1,n);
for i=1:n
    C{i}=f(i);
end

end

