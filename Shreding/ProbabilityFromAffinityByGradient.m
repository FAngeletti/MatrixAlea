function  [Pr,mu]= ProbabilityFromAffinityByGradient( Aff, widths, p ,Ws,mu,eps)
%ProbilityFromAffinityByGradient(Aff, widths, p ,Ws) Calculate the normalised probability from
%Affinity
%   Aff :positive nL*nPts matrice . nL number of shreds, nPts number of control points.
%   widths(k) : width of each section 
%   p(nPts) : probability function to be shred
%   Ws(nL) : weights of each shred function
%   The result Pr verify:
%       - sum (Pr(i,:) .*widths )=Ws(i)
%       - sum(Pr(:,k) ) = p(k)
%       - There is a vector mu such that
%            - Pr(i, k) = mu(i)* p(k) * Aff(i,k) / (sum(mu .* Aff(:,k) )

s=size(Aff);
nl=s(1);
Pr=zeros(s);
norm=zeros(s(2),1);

    function reach=optim(mu)
        reach=zeros(nl,1);
           for k=1:nl
                Pr(k,:)= mu(k).*Aff(k,:);
           end
           norm=sum(Pr);
           for l=1:nl
           reach(l)=sum(p.*widths.*Pr(l,:)./norm)-Ws(l);
           end
    end



mu=wishfullOptimisation4(@optim,100,1e-7,mu);

for i=1:nl
Pr(i,:)= p.*widths.*Pr(i,:)./norm;
end



end

