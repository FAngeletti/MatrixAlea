addpath /Users/pabry/MATLAB/UTILS_STAT/

%Dimension de la matrice
d=64;


%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d)./d;

v=1;
f=@(x) lnormal(x,v);
L0=lnormal(0,v);


%Matrice de lois
L1=cell(d,d);

% Cycle long
for i=1:d
    for j=1:d
		L1{i,j}=L0;
    end
end

for i=1:(d)
	L1{i,i}=f(i);
	L1{i,mod(i,d)+1}=f(i);
	L1{i,1+ mod( (i-1) +d -1 ,d)}=f(i);
end



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
Ts=CircularTs(lambdas)




% taille de la série temporelle
n=10000;

Law1=matrixLaw(A, E , L1 , n);

% Matrice de moyenne
M1=Law1.matMq(1);

% Coefficients associés à chaque valeurs propres
coeffs=CircularEigenCoeffs(M1)

len=400;





x1=Law1.rv();

figure(1); clf;
plot(x1);

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
