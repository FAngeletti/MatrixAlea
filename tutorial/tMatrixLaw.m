% Desired size of the signal
n=1000;

% Projection matrix
A = [ 1 1; 1 1];

% Structure matrix
p=0.8; q= 0.2;
E = [ p q; q p];

% Distribution matrix
Laws={ lNormal(-1,1), lNormal(0,2); lNormal(1,1), lNormal(0,2)};

% matrix Law
mL= matrixLaw(A,E,Laws,n);

% Signal synthesis
X=mL.rv();

plot(X)
