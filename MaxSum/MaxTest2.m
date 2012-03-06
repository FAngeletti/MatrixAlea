d=2;
p=1;
r=1;
E= [p r; 0 p];
A= ones(d);
L1=lgamma(1,1);
L2=lTranslated(lgamma(1,1),1);
L3=lnormal(1,1);

P= { L1 L1; L2 L2};

n=1000;
Lm=matrixLaw(A,E,P,n);


r=10000;

Vs=zeros(r,1);
for i=1:r
    x=Lm.rv();
 
    Vs(i)=max(Lm.rv())-log(r);
end


hist(Vs)
