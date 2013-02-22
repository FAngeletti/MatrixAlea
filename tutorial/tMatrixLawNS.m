% Desired size of the signal
n=1000;

% Projection matrix
A = [ 1 1; 1 1];

% Structure matrix
p=0.8; q= 0.2;
E = [ p q; q p];

% Distribution matrix
Laws{1} = { lNormal(-1,1), lNormal(0,2); lNormal(1,1), lNormal(0,2)};
Laws{2}=  {  lGamma(1,1), lGamma(2,1) ; lGamma(2,1),  lGamma(1,2)};

fLaws=@(n) Laws{n};

% matrix Law
mL= matrixLawNS(A,E,fLaws,n);

% Signal synthesis
X=mL.rv();

plot(X)
