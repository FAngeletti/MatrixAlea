function lambdas=CircularEigenValues(E)

P= E(1,:);
n=length(P);

nvp=n-1;


%Calcul des longueurs de corrélations et pseudo période
lambdas=zeros(1,nvp);

for k=1:nvp 
	lambdas(k)=evalPoly(P,k);
end


end

function w=w(d,k)
	i=complex(0,1);
	w= exp( i *(2 *pi*k) ./ d);
	
end

function y=evalPoly(P,k)
d=length(P);

tmp= w(d,k*(0:(d-1)) ) ;
y=sum(P.*tmp);
end
