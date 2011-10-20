

%Dimension de la matrice
d=2;


%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d);


% paramètre des gaussiennes;
%moy=zeros(d);
%var=ones(d);

moy= [0 0;
      0 0 ];

var= [ 1 3; 
       3 2 ];

%Matrice de lois
L=cell(d,d);

for i=1:d
for j=1:d
L{i,j}=lnormal(moy(i,j),var(i,j) ) ;
end
end

%Matrice de structure
p=0.98;
q=1-p;

E=[ p q;
    q p ];

% taille de la série temporelle
n=4000;

Law=matrixLaw(A, E , L , n);

x=Law.rv();

figure(1); clf;
plot(x);


[hh,bh,gh]=hist1d(x,100);
pdftheo=Law.pdf(bh);

figure(2) ; clf 
  plot(bh,hh); hold on ; plot(bh,gh,'--') ; grid on ; 
  plot(bh,pdftheo,'r--') ; 

xc = xcov(x,'coeff') ;
xc=xc(n:end);
ytheo= Law.corrVect();

len=100;
 
figure(3) ; clf 
plot(xc(1:len) ) ;  grid on; hold on;
plot(ytheo(1:len),'--r');

