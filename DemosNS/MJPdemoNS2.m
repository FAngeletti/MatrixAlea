%% Définition des lois sources
nP=2;
%Parents={lgamma(1,3), lgamma(1,3)} ;
Parents={lnormal(0,1), lnormal(0,1) };
figcount=1;

%% Découpage des sous lois
%Nombre de morceaux lois
nL=cellgen(@(i) 2, nP);


% Ws{i}(k) Poids associé à la k-ème sous-lois de la loi source Parents{i}
Ws=cellgen(@(i) ones(nL{i},1),nP );


%Définitions des points de contrôle par lois sources
% nps{i} nombre de points de contrôle du décupage de la loi parente Parents{i}
nps=cellgen( @(i) 1000,nP);

%Choix des points de contrôles
% pts{i} points de contrôle de la loi Parents{i}
% Choisi aléatoirement ici pour respecterapproximativement le découpage en
% quantile. ( les points auraient été de toutes façon triés, et les extrema de la
% loi parente insérés )
pts=cellgen(@(i) sort(vRandL(Parents{i},nps{i})), nP);
%pts=[-2,-1,0.5,0.75,1,2];

%Définition des fonctions d'affinité
% la fonction d'affinité funAff{i}{j}(x) correspond à l'affinité épprouvé par
% de la j-ème sous loi de la loi parente i au point x. L'attribution des
% densités de probabilités dans la fonction SchredByAffinity est faites en
% essayant de préserver les ratios d'affinités tout en respectant les
% contraintes de normalisation


%Nous allons les définir à partir d'un noyau et de points de contrôles

% Le noyau
kernel=@(v,x) exp(-(x-v(1)).^2./(2*v(2)^2));
%kernel=@(v,x) exp(-abs(((x-v(1))./v(2))));

%kernel=@(v,y) 1./(1+((v(1)-y)/v(2)).^2); 
%kernel=@(v,x) (sinc((x-v(1))./v(2))).^2 ;
%kernel= @(v,x) (1.5+cos(2*pi*(x-v(1))))./(1+ ((v(1)-x)./v(2)).^2 )  ;
%kernel = @(v,x) 1+tanh((x-v(1))/v(2));
% Les points de contrôles

%Expression fonctionnel
%Kpts=cellgen(@(i) pts{i}(1:((pts{i}(end)-pts{i}(1) )/(nL{i}-1) ) :  

%Expression déplié
Kpts=cell(1, nP);
%for i=1:nP
%    Kpts{i}=zeros(nL{i}, 2);
%    qs= (1:nP)./(nP+1);
%    stride= (nps{i}-1)/(nL{i}-1) ;
%    Kpts{i}(:,1) = pts{i} (fix( nps{i} *qs) );
%    Kpts{i}(:,2)= 0.25;
%end

Kpts{1} = [ 0 1; 0 10];
Kpts{2} = [0 1; 0 10];

% Les fonctions d'affinités {i}{j} sont alors directement définis en
% centrant les noyaux autour des points de contrôle kpts{i}{j}
%funAff=cellgen( @(i) AffinityFromKernel(kernel,Kpts{i}), nP );

f=@(x) exp(-0.25*(x.*x));
g=@(x) exp( -(x.^2-2.5).^2);
fg={f g};
funAff={fg fg};


% À partir de ces données ShredByAffinity construit les nL{i} sous lois
% associés aux nP lois parentes.
% Laws{i}{k} correspond à la k-ème sous-lois associés à la i-ème loi
% parente
Laws=cellgen( @(i) ShredByAffinity(Parents{i},pts{i}, Ws{i}, funAff{i}), nP) ;


%% Définition de E
d=2;

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

p=0.01;
q=1-p;

% Matrice de structure E
E=p*Id+q*Jd;

%% Définition de A
A=ones(d);

%% Définition des  matrices de proabilité
% Définition de la fonction probabilité matriciel
% MP(k) représente la matrice de probabilité au temps k.
MP= @(i) {Laws{i}{1} Laws{i}{2}; Laws{i}{1} Laws{i}{2} } ;

%% Définition de la loi matricielles
Law=matrixLawNS(A,E,MP,nP);

%% Paramètre de plot
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 

%% Synthesis

% s : nombre de réalisations indépendantes
s=10000;

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

vy=-3:0.15:3;
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
  