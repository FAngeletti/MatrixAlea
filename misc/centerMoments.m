function mc=centerMoments(m)
%centerMoments : compute the centered moments from a matrix of moments
% m(k,q). For now, the moments is supposed to be an integer

s=size(m);
mc=m;
nM=s(2);

for q=1:nM
	cnp=1;
	for k=1:(q-1)
		mc(:,q)= mc(:,q) + ((-1)^k) * cnp.*mc(:,k).* mc(:,1).^(q-k);
		cnp= cnp*(q-k+1)/ k;  
	end

end

end
