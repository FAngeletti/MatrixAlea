
% dimension of the matrices
d=2

% size of the random series
n=10000 ;

%Structure matrice
p=0.98;
q=1-p ;
E= [p q; q p]

% Projection matrix
A=[1 1 ; 1 1 ]


% Parent Law
pL=lLaplace(1/sqrt(2));

% Discretization points
pts=-5:0.01:5;

%Weights
Ws=[1;1]


%Kernel
kernel=@(v,x) exp(-abs((x-v(1))./v(2)));

%Targeted moments
tmoments=[-0.1 0.2; 0.1 1.8 ];

% Initial guess for the kernel parameter
kstart=tmoments

% Définition des lois sectionnées
Ls=ShredWithKernel(pL , tmoments, Ws , pts, kernel, kstart);


Lm=Ls{1};
Lp=Ls{2};



% Main diagonal of the probability matrix
diag={Lp,Lm};

%Probability matrix
Ls=cell(d,d);

for i=1:d
		Ls{i,i}=diag{i};
     		Ls{i,1+mod(i,d)}=diag{i};
end

Law=matrixLaw(A,E,Ls,n);


sx=Law.rv();

figure(1);
plot(sx);

figure(2);
x=-5:0.25:5;
[b, pos]=hist(sx,100);
b= b ./ ( n* ( pos(2)-pos(1) ) );
bar(pos,b); hold on; plot(x, pL.pdf(x), 'r--'); hold off;


