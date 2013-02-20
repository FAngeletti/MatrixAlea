
%Dimension de la matrice
d=32;


%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d)./d;

%Matrice de structure
p=0.95;
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
targetCoeffs=zeros(1,fix(d/2));
targetCoeffs(1)=0.05;
targetCoeffs(3)=0.025;
targetCoeffs(5)=0.1;
targetCoeffs(8)=0.025;

if mod(d,2)==0
symCoeffs= [ 0 targetCoeffs fliplr(targetCoeffs(1:(end-1)))  ];
else
symCoeffs= [ 0 targetCoeffs fliplr(targetCoeffs)  ];

end

%Calcul des coefficients de la matrice
diag= fft(sqrt(symCoeffs));
%normalisation
diag=real(diag-sum(diag));


L0=lnormal(0,1);

lvar=0.8;
kernel= @(v,x) exp(-((x-v(1))./lvar).^2.);
%kernel= @(v,x) exp(-abs((x-v(1)))./lvar );
%kernel= @(v,x) exp(-((x-v(1))./var).^4.);
%kernel=@(v,x) 1./(1+ ((x-v(1)).^2)./lvar);

kstart=diag' ;

Ws=ones(length(diag),1); 
if(pts==0)
pts=DiscreteApproximationBase(L0,0.000001,-5,5);
end

Laws=ShredWithKernel(L0, diag', Ws,pts,kernel,kstart);

TestMixture(L0,Laws,Ws,10,50);
pdflen=max(abs(diag))/p + 4 ;


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
 %   L1{i,1+mod(i,d)}=Laws{i};
end



% taille de la série temporelle
n=5000;

Law1=matrixLaw(A, E , L1 , n);

% Matrice de moyenne
M1=Law1.matMq(1);

% Coefficients associés à chaque valeurs propres, calculés numériquement
coeffs=CircularEigenCoeffs(M1)

x1=Law1.rv();

figure(4); clf;
plot(x1);

figure(5);clf
xpos= (-pdflen):0.1:pdflen;
plot(xpos,Law1.pdf(xpos));

var= Law1.variance();

xc = xcov(x1,'coeff') ;
xc=xc((n+1):end);

xctheor=0*x1;
for t=1:length(x1)
	xctheor(t)=CircularCov(lambdas,coeffs,t);
end

xctheor=xctheor./var;
 
figure(6) ; clf 
plot(xc ) ;  grid on; hold on;
plot(xctheor, 'r--');
