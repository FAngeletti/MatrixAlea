function C= cellgen( f, sh )
%cellgen(f,sh) Create a cell C of shape sh, where C{i}=f(i); 


dim=length(sh);
if(dim==1)
    C=cell(dim,1);
else
    C=cell(sh);
end

indice=init(dim);

len =prod(sh);
for i=1:len
    C{indice{:}}=f(indice{:});
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

