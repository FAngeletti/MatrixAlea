addpath /Users/pabry/MATLAB/UTILS_STAT/

%Dimension de la matrice
d=32;


%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d);


f=@(x) lnormal(x,2);
L0=lnormal(0,1);

%Matrice de lois
L1=cell(d,d);

% Cycle long
for i=1:d
    for j=1:d
	if(i==j)
		L1{i,j}=f(i);
	elseif (abs(i-j)==1)
     		L1{i,j}=f(i) ;
	else
		L1{i,j}=L0;
	end
    end
end

%Matrice de structure
p=0.75;
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
%E= p.*Id+q.*(J) ;

%Structure symétrique
E=p.*Id+0.5*q*(J+J^(d-1) );
% taille de la série temporelle
n=10000;

Law1=matrixLaw(A, E , L1 , n);

x1=Law1.rv();

figure(1); clf;
plot(x1);



[hh,bh,gh1]=hist1d(x1,100);
pdftheo=Law1.pdf(bh);

figure(3) ; clf 
  plot(bh,hh,'k'); hold on ; grid on ; 
  plot(bh,pdftheo,'--') ; 
  


xc = xcov(x1,'coeff') ;
xc=xc((n+1):end);


len=1000;
ytheo1=zeros(len,1);

for i=1:len
ytheo1(i)=Law1.corr(i+1);
end
 
figure(4) ; clf 
plot(xc(1:len) ) ;  grid on; hold on;
plot(ytheo1,'r--');

