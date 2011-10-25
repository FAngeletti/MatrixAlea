%% Définition des lois sources
nP=2;
%Parents={lgamma(1,3), lgamma(1,3)} ;
Parents={lnormal(0,1), lgamma(2,1) };
figcount=1;

%% Découpage des sous lois
%Nombre de morceaux lois
nL={8 8};


% Ws{i}(k) Poids associé à la k-ème sous-lois de la loi source Parents{i}
Ws=cellgen(@(i) ones(nL{i},1),nP );



%Choix des points de contrôles
pts={(-4:0.01:4) (0:0.01:10) };
%On récupère le nombre de points de contrôle
npts=cellgen(@(i) length(pts{i}), 2 );

%On définit un noyau de découpe de loi
kernel=@(v,x) exp(-(x-v(1)).^2.);

% On choisit des moments cibles
Tmoments= { [-1.4;-0.7;-0.5;-0.2;0.2;0.5;0.7;1.4], [0.7;1.;1.3;1.6;1.9;2.5;3.;4] };
% Tmoments{k}(l,q) doit correspondre au q-ème moment de la sous-loi l de la loi parente k.
% Attention! La somme des moments partiels doit rester égale au moment de la loi parente.
% De plus, les bornes présentés, ci-haut sont très proches des bornes maximales. Les dépasser empêchera
% probablement Shred...Wishfully de converger. 


%En première approximation, on considère que les centres des noyaux devraient être proche des moments choisi
Kstart=Tmoments;


% À partir de ces données ShredWithKernelWishfully construit les nL{i} sous lois
% associés aux nP lois parentes.
% Laws{i}{k} correspond à la k-ème sous-lois associés à la i-ème loi
% parente

Laws=cellgen( @(i) ShredWithKernelWishfully(Parents{i} , Tmoments{i}, Ws{i}, pts{i}, kernel, Kstart{i} ), nP) ;


%% Définition de E
d=8;

% Matrice identité
Id=zeros(d);
for i=1:d
    Id(i,i)=1;
end

% Matrice circulaire
Jd=zeros(d);
for i=1:d
    Jd(i, mod(i, d) +1)=1;
end

% Matrice antidiagonale
Ad=zeros(d);
for i=1:d
    Ad(i, d-i+1)=1;
end

p=0.01;
q=1-p;

% Matrice de structure E
%E=p*Id+q*Jd;
E=p*Id+q*Ad;



%% Définition de A
A=ones(d);

%% Définition des  matrices de proabilité
% Définition de la fonction probabilité matriciel
% MP(k) représente la matrice de probabilité au temps k.
Mc=cell(nP,1);

for i=1:nP
Mc{i}=cell(d,d);

for k=1:d
    for l=1:d
        Mc{i}{k,l}=Parents{i};
    end
end

for k=1:d
    Mc{i}{k,k}=Laws{i}{k};
    Mc{i}{k,d-k+1}=Laws{i}{k};
end
end

MP= @(i) Mc{i} ;

%% Définition de la loi matricielles
Law=matrixLawNS(A,E,MP,nP);

%% Paramètre de plot
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 

%% Synthesis

% s : nombre de réalisations indépendantes
s=1000;

%Synthèse de s réalisations de Law
x=zeros(nP,s);
for k=1:s
    tmp=Law.rv();
   x(:,k)=tmp;
end

% Min et max
 xmin=min(x(1,:)); xmax=max(x(1,:));
 ymin=min(x(2,:)); ymax=max(x(2,:));

% Affichage des réalisations composantes par composantes
for i=1:nP
figure(figcount); clf;
plot(x(i,:),'k');grid on  
set(gca,'FontSize',fontsize) 
h(1) = ylabel('Data') ; 
h(2) = title(sprintf('X_%d',i)) ; 
set(h,'FontSize',fontsize) ; 
figcount=figcount+1;
end

figure(figcount); clf;
plot(x(1,:), x(2, :) ,'.');
figcount=figcount+1;

%%  pdf uni

% Vérification des pdf univariés
for i=1:nP
    
[hh,bh,gh1]=hist1d(x(i,:),50);
pdftheo=Law.pdf(bh,i);
pdfTarget=Parents{i}.pdf(bh);

  figure(figcount) ; clf 
  plot(bh,pdftheo(:,1),'b--','LineWidth',linewidth,'MarkerSize',markersize);grid on; hold on;
  plot(bh,pdfTarget,'r--','LineWidth',linewidth,'MarkerSize',markersize);
plot(bh,hh,'k','LineWidth',linewidth,'MarkerSize',markersize ) ;  

color=hsv(nL{i});
for k=1:nL{i}
    plot(bh,Laws{i}{k}.pdf(bh)./(nL{i}),'color',color(k,:),'Linestyle', ':' ,'LineWidth',linewidth,'MarkerSize',markersize);
end

set(gca,'FontSize',fontsize) 
axis([xmin xmax 0 0.55])
h(1) =  ylabel('Marginal') ;  
set(h,'FontSize',fontsize) ;
figcount=figcount+1;
end


%% pdf multivarié theorique

nx=50;
ny=50;
 
%stepx= (xmax - xmin) /(nx-1);
%stepy= (ymax - ymin) /(ny-1);
 

%vy=ymin:stepy:ymax;
%vx=xmin:stepx:xmax;

vy=0:0.25:8;
vx=-3:0.15:3;

nx=length(vx);
ny=length(vy);


 [Gx,Gy]=meshgrid(vx,vy);
 
 Z=zeros(ny,nx);
 

 
 for i=1:nx
     for j=1:ny
        Z(j,i)=(Law.mvPdf([vx(i) vy(j)])) ;
     end
 end
figure(figcount); clf;
 surf(Gx,Gy,Z);
 
  %% Corr
  covTheor=Law.mvMoments([1 2], [1 1]) - Law.moments(1,1)*Law.moments(1,2)
  corrTheo=covTheor/sqrt( (Law.moments(2,1)-Law.moments(1,1)^2)*(Law.moments(2,2)-Law.moments(1,2)^2) )
  
  cov=sum(x(1,:).*x(2,:))/s - ( sum(x(1,:))*sum(x(2,:))/(s^2) )
  covCorrcoef=corrcoef(x')
  
  %% Square Correlation
  x2=x.^2;
  cov2Theor=Law.mvMoments([1 2], [2 2]) - Law.moments(2,1)*Law.moments(2,2)
  corr2Theo=cov2Theor/sqrt( (Law.moments(4,1)-Law.moments(2,1)^2)*(Law.moments(4,2)-Law.moments(2,2)^2) )
  
  cov2=sum(x2(1,:).*x2(2,:))/s - ( sum(x2(1,:))*sum(x2(2,:))/(s^2) )
  corr2Corrcoef=corrcoef(x2')
  
