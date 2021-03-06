% Desired size of the signal
n=1000;

% Projection matrix : Initial distribution of the hidden Markov chain
A = [ 1 1; 1 1];

% Structure matrix : Transition matrix of the hidden Markov chain
p=0.8; q= 0.2;
E = [ p q; q p];

% Distribution matrix: Matrix of the laws of the observable
sigma=1;
Laws={ lNormal(-1,sigma), lNormal(0,2*sigma) ; 
       lNormal(1,sigma),  lNormal(0,2*sigma) };

% matrix Law
mL= matrixLaw(A,E,Laws,n);

% Signal synthesis
X =mL.rv();

plot(X)
