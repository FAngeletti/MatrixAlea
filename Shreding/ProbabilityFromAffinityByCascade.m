function  Pr= ProbabilityFromAffinityByCascade( Aff, widths, p ,Ws)
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
Pr=zeros(s);
disp(Ws)
for it=1:(nl-1)
    P=sum(p.*widths); %disp(P);
    AffStar= sum(Aff((it+1):end,1:end),1);
    AffSel=Aff(it,:);
    target=P-Ws(it); %fprintf(1,'target=%f\n',target); 
    
    f=@(mu) sum( (p.*AffStar.*widths) ./ (mu * AffSel + AffStar) ) -target;
    mu=fzero(f,[0 100]); fprintf(1,'Step %d : mu=%f\n', it, mu);
    
    tmp=(mu*AffSel+ AffStar);
    Pr(it,:)= (mu*AffSel.* p) ./tmp; 
   
    p=p-Pr(it,:);
end

Pr(nl,:)=p;


end

