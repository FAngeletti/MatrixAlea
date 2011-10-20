addpath /Users/pabry/MATLAB/UTILS_STAT/

%Dimension de la matrice
d=8;


%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d)./d;

%Matrice de structure
p=0.9;
q=1-p ;

J=zeros(d);
for(i=1:(d-1))
J(i,i+1 )=1;
end 
J(d,1)=1;

Id=zeros(d);
for(i=1:d)
Id(i,i)=1;
end 


% Structure asymétrique
E= p.*Id+q.*(J) ;

%Structure symétrique
%E=p.*Id+0.5*q*(J+J^(d-1));

% Valeurs propres associées à la matrice de structure
lambdas= CircularEigenValues(E)

% Longueurs de corrélation et pseudo-périodes associées 
taus= CircularTaus(lambdas)

corrlen=max(taus);
len=round(5*corrlen);

Ts=CircularTs(lambdas)

% Coefficients associés aix lambda_j souhaités
targetCoeffs= [ 1.5 0 0 0];

if mod(d,2)==0
symCoeffs= [ 0 targetCoeffs fliplr(targetCoeffs(1:(end-1)))  ];
else
symCoeffs= [ 0 targetCoeffs fliplr(targetCoeffs)  ];

end

%Calcul des coefficients de la matrice
diag= fft(sqrt(symCoeffs));
%normalisation
diag=diag-sum(diag)

var=2;
L0=lnormal(0,var);

%kernel= @(v,x) exp(-((x-v(1))./var).^2./2);
%kernel= @(v,x) exp(-2*abs((x-v(1))) );
kernel= @(v,x) exp(-((x-v(1))./var).^4.);


Ws=ones(length(diag),1); pts=[-8:0.025:8];
Laws=ShredWithKernelWishfully(L0, diag', Ws,pts,kernel,diag');

TestMixture(L0,Laws,Ws,1000,50);
pdflen=max(abs(diag))/p + 4 * var;


%Matrice de lois
L1=cell(d,d);

% Cycle long
for i=1:d
    for j=1:d
		L1{i,j}=L0;
    end
end

for i=1:(d)
	L1{i,i}=Laws{i};
end



% taille de la série temporelle
n=10000;

Law1=matrixLaw(A, E , L1 , n);

% Matrice de moyenne
M1=Law1.matMq(1);

% Coefficients associés à chaque valeurs propres, calculés numériquement
coeffs=CircularEigenCoeffs(M1)

x1=Law1.rv();

figure(4); clf;
plot(x1(1:(5*len)));

figure(5);clf
xpos= (-pdflen):0.1:pdflen;
plot(xpos,Law1.pdf(xpos));

var= Law1.variance();

xc = xcov(x1,'coeff') ;
xc=xc((n+1):end);

xctheor=zeros(1,len);
for t=1:len
	xctheor(t)=CircularCov(lambdas,coeffs,t);
end

xctheor=xctheor./var;
 
figure(6) ; clf 
plot(xc(1:len) ) ;  grid on; hold on;
plot(xctheor, 'r--');
