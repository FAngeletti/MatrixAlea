function  [Pr,mu]= ProbabilityFromAffinityByCascade( Aff, widths, p ,Ws,mu,eps)
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
eps=0.00001;
for it=1:(nl-1)
    P=sum(p.*widths); %disp(P);
    AffStar= sum(Aff((it+1):end,1:end),1);    
    AffSel=Aff(it,:);
    target=P-Ws(it); %fprintf(1,'target=%f\n',target); 
    f=@(x) sum( (p.*AffStar.*widths) ./ (x .* AffSel + AffStar) ) -target;


    bsup=(P/target-1)*max(AffStar./AffSel);
    
    mu(it)=fzero(f,[eps bsup]); %fprintf(1,'Step %d : mu=%f\n', it, mu(it));
    
    tmp=(mu(it)*AffSel+ AffStar);
    Pr(it,:)= (mu(it)*AffSel.* p) ./tmp; 
   
    p=p-Pr(it,:);
end

Pr(nl,:)=p;


end

