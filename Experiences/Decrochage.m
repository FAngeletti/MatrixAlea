%Loi parente
L0=lnormal(0,1)

%dimension
d=4;

%% Découpage des lois
kernel=@(v,x) exp(-(x-v(1)).^2);
pts=-4:0.01:4;
Ws=ones(d,1);
tmoments=[-1.2; -0.3; 0.3 ; 1.2 ];
kstart=tmoments;

Laws=ShredWithKernelWishfully(L0, tmoments,Ws,pts, kernel, kstart);
%Test du découpage
%TestMixture(L,Laws,Ws,1000,50);

%% Definition de MP

MP=cellgen(@(i,j) L0, [d d]);

%% Choix de p
for i=1:d
    MP{i,i}=Laws{i};
end

% vecteurs de moments
Am=arraygen( @(i) MP{i,i}.moments(1), d );
Vm=arraygen( @(i) MP{i,i}.moments(2), d );

kappa= sum(Vm)./sum(Am.^2)

r=(1-1/kappa)^(1/d)
p=1-1/r
q=1-p

%% Definition de A
A=ones(d);

%% Définition de E

% Matrice identité
Id=zeros(d);
for i=1:d
    Id(i,i)=1;
end

% Matrice circulaire
Jd=zeros(d);
for i=1:d
    Jd(i, mod(i, d) +1)=1;
end

E=p*Id+q*Jd;

n=1000;
L=matrixLawNS(A,E,MP,n);

x=L.rv();

figure(1);
cov=xcov(x,'coeff');
plot(cov,'k'); hold on; frid on;
av2=L.moments(1)^2;
var=L.moments(2) -av2;
covt=arraygen( @(i) L.mvMoments( [1 i], [1 1] )-av2);
covt= [1 covt./var];
plot(covt,'r-'); hold off

