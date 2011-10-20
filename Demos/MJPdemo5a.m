addpath /Users/pabry/MATLAB/UTILS_STAT/

%Dimension de la matrice
d=6;


%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d);


% paramètre des gaussiennes;
%moy=zeros(d);
%var=ones(d);

L0=lnormal(0,1);
Lp=lnormal(1,2);
Lm=lnormal(-1,2);

% Même marginale
% Différence de corrélation
diag1={Lp Lp Lp Lm Lm Lm};
diag2={Lp Lm Lp Lm Lp Lm};

%Matrice de lois
L1=cell(d,d);
L2=cell(d,d);


for i=1:d
    for j=1:d
	if(i==j)
		L1{i,j}=diag1{i};
		L2{i,j}=diag2{i};
	else
     		L1{i,j}=L0 ;
		L2{i,j}=L0;
	end
    end
end

%Matrice de structure
p=0.98;
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

E= p.*Id+q.*(J) ;

% taille de la série temporelle
n=20000;

Law1=matrixLaw(A, E , L1 , n);
Law2=matrixLaw(A,E,L2,n);

x1=Law1.rv();
x2=Law2.rv();

figure(1); clf;
plot(x1);

figure(2); clf;
plot(x2);






[hh,bh,gh1]=hist1d(x1,100);
pdftheo=Law1.pdf(bh);

figure(3) ; clf 
  plot(bh,hh,'k'); hold on ; grid on ; 
  plot(bh,pdftheo,'--') ; 
  
[hh,bh,gh1]=hist1d(x2,100);
pdftheo=Law2.pdf(bh);

 plot(bh,hh,'b');
 plot(bh,pdftheo,'g--') ; 


xc = xcov(x1,'coeff') ;
xc=xc((n+1):end);

xc2 = xcov(x2,'coeff') ;
xc2=xc2((n+1):end);

len=300;
ytheo1=zeros(len,1);
ytheo2=zeros(len,1);

for i=1:len
ytheo1(i)=Law1.corr(i+1);
ytheo2(i)=Law2.corr(i+1);
end
 
figure(4) ; clf 
plot(xc(1:len) ) ;  grid on; hold on;
plot(ytheo1,'r--');
plot(xc2(1:len),'k' ) ;  
plot(ytheo2,'g--');
