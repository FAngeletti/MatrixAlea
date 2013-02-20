function Laws = ShredWithKernel(L , moments, Ws0, pts0, Kern, Kstart )
%ShredWithKernelWishfully(L , momentsTarget, Ws, pts0, Kern, Kstart ) Try to
%shred the law L into shred Laws of weigths Ws such that E Laws{i} ^ q = moments(i,q),
%using control point pts and the kernel Kern starting with kernel
%parameters Kstart.

Ws=Ws0./sum(Ws0);

s=size(Kstart);
nL=s(1);
dimK=s(2);
sM=size(moments);
dimM=sM(2);

if( sM(1) ~= nL )
	error('MJP:IncoherentWish',  'Required moments for %d laws, however the kernel initial parameters provides parameter for %d laws',sM(1), nL) ;
end

if(dimK<dimM)
    warning('Shred:DimWarning','Trying to fit kernel shred to moments with less parameter than kernel');
end

for q=1:dimM
    Mt=sum(Ws.*moments(:,q));
    if( abs(Mt-L.moments(q)) >nL*eps )
        error('MJP:ImpossibleWish', 'Required moments (sum : %d) are incoherent with the parent law (%f)',Mt,L.moments(q) );
    end
end

%Construction of the remainder law
R=lRemainder(L,pts0);
Rq=zeros(1,dimM);


Wr=R.weight();

pts=R.points();
p=R.dps;
larg=diff(pts);
center=(pts(2:end)+pts(1:(end-1) ))./2;
npts=length(pts);

pR=ChooseMixedShredLaw(R.weight(), Ws);
for q=1:dimM
Rq(q)=Wr*R.moments(q);
end

Ws(pR)=Ws(pR)-Wr;



Aff=zeros(nL,npts-1);
Cmom= zeros(nL,dimM);
Pr=0.*Aff;
precompQ=zeros(dimM,npts-1);
for q=1:dimM
precompQ(q,:)= diff(pts.^(1+q))./(q+1);
end
cmoments=centerMoments(moments);
mu=ones(nL,1);
intEps=1e-9;
% Definition of the function to be optimised
    function fitness=Optim(Params)
%Construction of the affinity matrice
        for j=1:nL
            Aff(j,:) = Kern(Params(j,:), center);
        end
	
        
    %[Pr,mu]=ProbabilityFromAffinityByDyadicDivision(Aff, larg, p, Ws,mu,intEps) ;
    %  [Pr,mu]=ProbabilityFromAffinityByCascade(Aff, larg, p, Ws,mu,intEps) ;
     [Pr,mu]=ProbabilityFromAffinity(Aff,larg,p,Ws,mu,intEps);
  

      for qi=1:dimM
          for k=1:nL
              Cmom(k,qi) = sum(Pr(k,:).* precompQ(qi,:))./Ws(k);
          end
      end
      
      Cmom(pR,:) =(Ws(pR).*Cmom(pR,:)+ Rq)./(Ws(pR)+Wr);  
      fitness=centerMoments(Cmom)-cmoments;
     % fitness=reshape(Cmom,dimM*nL,1);
%Computation of the associated moments      
    end

param=Kstart;
epsS1=1e-3;
epsT=1e-6;
options=optimset('Tolfun',epsT, 'TolX',1e-8);

[param,err]=AdhocOptimisation4(@Optim,100, epsT, param);
%[param,err]=fsolve(@(x) reshape(Optim(x),dimM*nL,1), param,options);

disp('param:');
 disp(param);
 disp('error:');
 disp(err);
   





   
  
  
 
    
    Laws=cell(1,nL);
%Constructing the laws.
for i=1:nL
    L=lSimpleF(pts,Pr(i,:));
    if (i==pR)
        Laws{i}= lMixture({L R}, [Ws(i) Wr] );
    else
        Laws{i}=L;
    end
end

end

