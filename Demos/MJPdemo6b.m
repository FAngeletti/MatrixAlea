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
targetCoeffs= [ 0.5 0 2 1];

if mod(d,2)==0
symCoeffs= [ 0 targetCoeffs fliplr(targetCoeffs(1:(end-1)))  ];
else
symCoeffs= [ 0 targetCoeffs fliplr(targetCoeffs)  ];

end


diag= fft(sqrt(symCoeffs));


%for k=1:(d/2)
%		diag=diag+ sqrt(2*targetCoeffs(k)) * cos( (tic *k ) * (1:d) );
%end


% Afin d'obtenir la corrélation souhaité, il est nécessaire que <f(x)>=x
v=2;
f=@(x) lnormal(x,v);
L0=lnormal(0,v);
pdflen=max(abs(diag))/p + 4 * v;


%Matrice de lois
L1=cell(d,d);

% Cycle long
for i=1:d
    for j=1:d
		L1{i,j}=L0;
    end
end

for i=1:(d)
	m=f(diag(i)/p );
	L1{i,i}=m;
%	L1{i,mod(i,d)+1}=m;
%	L1{i,1+ mod( (i-1) +d -1 ,d)}=m;
end



% taille de la série temporelle
n=10000;

Law1=matrixLaw(A, E , L1 , n);

% Matrice de moyenne
M1=Law1.matMq(1);

% Coefficients associés à chaque valeurs propres, calculés numériquement
coeffs=CircularEigenCoeffs(M1)

x1=Law1.rv();

figure(1); clf;
plot(x1(1:(5*len)));

figure(2);clf
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
 
figure(4) ; clf 
plot(xc(1:len) ) ;  grid on; hold on;
plot(xctheor, 'r--');
