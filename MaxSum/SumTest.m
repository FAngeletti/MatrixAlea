d=4;
p=1;
r=1;
E= [p r/2 r/2 0; 0 p 0 r; 0 0 p r; 0 0 0 p];
A= ones(d);
L1=lnormal(0,1);
L2=lnormal(-1,1);
L3=lnormal(1,1);

P= { L1 L1 L1 L1; L2 L2 L2 L2; L3 L3 L3 L3; L1 L1 L1 L1};

n=1000;
Lm=matrixLaw(A,E,P,n);


r=1;

Vs=zeros(r,1);
for i=1:r
    x=Lm.rv();
    plot(x);
    Vs(i)=sum(Lm.rv());
end

