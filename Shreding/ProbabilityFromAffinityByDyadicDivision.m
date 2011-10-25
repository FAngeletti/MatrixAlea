function  [Pr,mu]= ProbabilityFromAffinityByDyadicDivision( Aff, widths, p ,Ws,mu)
%ProbilityFromAffinityByDyadic(Aff, widths, p ,Ws) Calculate the normalised probability from
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

if nl==1
    Pr=p;
else
pivot=fix(nl/2);
    
    iinf=(1:pivot);
    isup=((pivot+1):nl);
    P=sum(p.*widths); %disp(P);
    AffInf=Aff(iinf,:);
    AffSup=Aff(isup,:);
    Wsup=Ws(isup);
    Winf=Ws(iinf);
    
    AffSel= sum(AffInf,1);    
    AffStar=sum(AffSup,1);
    W=sum(Winf);
    
    
    target=P-W; %fprintf(1,'target=%f\n',target);
    f=@(x) sum( (p.*AffStar.*widths) ./ (x .* AffSel + AffStar) ) -target;
    
    bsup=(P/target-1)*max(AffStar./AffSel);
    mu=fzero(f,[0 bsup]); %fprintf(1,'Step %d : mu=%f\n', it, mu(it));
    
    tmp=(mu*AffSel+ AffStar);
    Pinf= (mu*AffSel.* p) ./tmp;
    Psup=p-Pinf;
    
    [Pr(iinf,:),mu]=ProbabilityFromAffinityByDyadicDivision(AffInf,widths,Pinf,Winf,mu);
    [Pr(isup,:),mu]=ProbabilityFromAffinityByDyadicDivision(AffSup,widths,Psup,Wsup,mu);

end



end

