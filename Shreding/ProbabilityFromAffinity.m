function  [Pr,mu]= ProbabilityFromAffinity( Aff, widths, p ,Ws,mu,eps)
%ProbilityFromAffinity(Aff, widths, p ,Ws) Calculate the normalised probability from
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
npts=s(2);
Pr=zeros(s);
%disp(Ws)
if( nargin<5)
mu=ones(nl,1);
end
pw=p.*widths;
affw=zeros(1,npts);
wm=zeros(nl,1);

    function er=errorfun(mu)
        
        for k=1:npts
        affw(k)= sum(abs(mu).* Aff(:,k) );
        end
        
        for i=1:nl
        wm(i)= sum( ( abs(mu(i)).*Aff(i,:).*pw )./ affw );
        end
        er=wm-Ws;
    end

options=optimset('Display','off','TolFun', eps);
mu=abs(fsolve(@errorfun, mu,options));

%mu=wishfulOptimisation4(@errorfun,100,1e-7,mu,0);

for j=1:nl
Pr(j,:)=mu(j).*p.*Aff(j,:) ./ affw;
end

end

