function [ v ] = vRandL(L,n )

v=zeros(1,n);
for i=1:n
    v(i)=L.rv();
end

end

