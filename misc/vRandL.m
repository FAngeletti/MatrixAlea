function [ v ] = vRandL(L,n )
%vRandL(L,n) Generate a vector of size n with law L  
v=zeros(1,n);
for i=1:n
    v(i)=L.rv();
end

end

