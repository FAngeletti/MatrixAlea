% Définition des lois par coupure

clear all
close all

%% Construction

% taille de la série temporelle
n=2 ;

%Matrice de structure
p=0.98;
q=1-p ;

% Lois de bases
Lg={lgamma(1,2),lgamma(2,3)};

% Cut
cut=cell(2,1);
for i=1:n
cut{i}= fzero( @(x) Lg{i}.cumulative(x)-0.5, [0,10]);
end

%Découpage des lois de bases
Lp=cell(2,1);
Lm=cell(2,1);
for i=1:n
Lp{i}=lTruncated(Lg{i},[0,cut{i}]);
Lm{i}=lTruncated(Lg{i},[cut{i},Inf]);
end

Ls=cell(n,1);

%Construction de la matrice
for i=1:n
    tmp={Lm{i} Lp{i}; Lm{i} Lp{i} };
    Ls{i}=tmp;
end

%Dimension de la matrice
d=2;

% affichage
len = 300 ; 
V = [0 len -0.1 1] ; 
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 

%Définition de l'opérateur de projection L(M) = <A,M> 
A=ones(d);

J=zeros(d);
for(i=1:(d-1))
J(i,i+1 )=1;
end 
J(d,1)=1;

Id=zeros(d);
for i=1:d
Id(i,i)=1;
end 

%matrice de structure 
E=p*Id+q*J;



gen= @(x)( @(n) x{n});

Law=matrixLawNS(A,E,gen(Ls),n);


%% Synthesis

s=1000;
x=zeros(n,s);



for k=1:s
    tmp=Law.rv();
   x(:,k)=tmp;
end
figcount=1;
for i=1:n
figure(figcount); clf;
plot(x(i,:),'k');grid on  
set(gca,'FontSize',fontsize) 
h(1) = ylabel('Data') ; 
h(2) = title(sprintf('X_%d',i)) ; 
set(h,'FontSize',fontsize) ; 
figcount=figcount+1;
end


%%  pdf uni


for i=1:n
    
[hh,bh,gh1]=hist1d(x(i,:),50);
pdftheo=Law.pdf(bh,i);
pdfTarget=Lg{i}.pdf(bh);
% figure(3) ; clf 
%   plot(bh,hh,'k'); hold on ; grid on ; 
%   plot(bh,pdftheo,'--') ; 
  
  figure(figcount) ; clf 
  plot(bh,pdftheo(:,1),'b--','LineWidth',linewidth,'MarkerSize',markersize);grid on; hold on;
  plot(bh,pdfTarget,'r--','LineWidth',linewidth,'MarkerSize',markersize);
plot(bh,hh,'k','LineWidth',linewidth,'MarkerSize',markersize ) ;  
set(gca,'FontSize',fontsize) 
axis([0 10 0 0.55])
h(1) = ylabel('Marginal') ; 
set(h,'FontSize',fontsize) ;
figcount=figcount+1;
end

 %% pdf multivarié theorique
v=0:0.2:10;
ngrid=length(v);
 [Gx,Gy]=meshgrid(v);
 
 Z=0.*Gx;
 
 for i=1:ngrid
     for j=1:ngrid
        Z(i,j)=(Law.mvPdf([v(i) v(j)])) ;
     end
 end
figure(5);
 surf(Gx,Gy,Z);
 
  %% Corr
  covTheor=Law.mvMoments([1 2], [1 1]) - Law.moments(1,1)*Law.moments(1,2)
  corrTheo=covTheor/sqrt( Law.moments(2,1)*Law.moments(2,2) )
  
  cov=sum(x(1,:).*x(2,:))/s - ( sum(x(1,:))*sum(x(2,:))/(s^2) )

