function A= arraygen( f, sh )
%arraygen(f,sh) Create an array A of shape sh, where C(i)=f(i); 

dim=length(sh);
if(dim==1)
    A=zeros(dim,1);
else
A=zeros(sh);
end

indice=init(dim);

len =prod(sh);
for i=1:len
    A(indice{:})=f(indice{:});
    indice=incr(sh,indice);
end

end

function c=init(dim)

c=cell(1,dim);

for i=1:dim
    c{i}=1;
end

end

function v=incr(sh,indice)
    v=indice;
    v{1}=v{1}+1;
    i=1;
    dim=length(sh);
    while(i < dim && v{i}>sh(i) )
        v{i}=1; v{i+1}=v{i+1}+1; i=i+1;
    end

    
end

