function [ Laws ] = ShredByAffinity( Ls, pts0 , Ws0 , funAff, verbose )
%ShredByAffinity(Ls, pts, Ws, funAff, verbose) Shred a parent distribution Ls in
%different partial law, in such a way that the mixture with weight Ws recompose the original law 
% Ls : Parent Law
% pts : Control point of the shreding
% Ws : wieghts of each partial law
% funAff(k) : affinity function of the partial law k.
if(nargin<5)
    verbose=1;
end
%Debugging function
    function dprintf(format,varargin)
        if(verbose==1)
            fprintf(1,format,varargin{:});
        end
    end

    
    function ddisp(x)
        if(verbose==1)
            disp(x);
        end
    end

% Determining the remainder law
R=lRemainder(Ls, pts0);

% Warning! The division algorithm needs to insert the extrema points of pdf
% amongst the points. The constructor of the remainder law construct this extended points, and we use it.  
pts=R.points();
np=length(pts);
Ws=Ws0./sum(Ws0);
nl=length(Ws);
larg=diff(pts);

center=(pts(2:end)+pts(1:(end-1)))/2; 

% Determining the affinity matrice
Aff=zeros(nl,np-1);
for k=1:nl
Aff(k,:)=funAff{k}(center);
end

%Choosing which law will be mixed with the remainder law
Wr=R.weight();dprintf('Remainder weigth : %f \n',Wr);
pR=ChooseMixedShredLaw(Wr,Ws);


Ws(pR)=Ws(pR)-Wr;

ddisp('Ws='); ddisp(Ws);

p=R.dps;
P=sum(p.*larg);
dprintf('probabilité totale disponible : %f \n', P);

% Solving the normalisation problems
Pr=ProbabilityFromAffinity(Aff,larg,p,Ws);

%Aff
if (verbose==1)
for k=1:nl
    fprintf(1,' Sum(Aff(%d,:) * larg) = %f \n', k , sum(Pr(k,:).*larg));
end
end

Laws=cell(nl,1);

%Constructing the laws.
for i=1:nl
    L=lSimpleF(pts,Pr(i,:));
    if (i==pR)
        Laws{i}= lMixture({L R}, [Ws(i) Wr] );
    else
        Laws{i}=L;
    end
end

end

