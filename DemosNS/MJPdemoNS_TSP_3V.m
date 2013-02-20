
figcount=1;
printf=@(s) @() print('-depsc',s);
printing=@(s) print('-depsc',s); 


fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 3 ; 

pplot=@(x,y,s) plot(x,y,s,'LineWidth',linewidth,'MarkerSize',markersize);

%% Définition des lois sources
nP=3;
%Parents={lgamma(1,3), lgamma(1,3)} ;
Parents={lnormal(0,1), lgamma(2,1), lgamma(1,2) };
figcount=1;

%% Découpage des sous lois
%Nombre de morceaux lois
nL={2 2 2};


% Ws{i}(k) Poids associé à la k-ème sous-lois de la loi source Parents{i}
Ws=cellgen(@(i) ones(nL{i},1),nP );



%Choix des points de contrôles
pts={(-4:0.01:4) (0:0.01:10) (0:0.01:10) };
%On récupère le nombre de points de contrôle
npts=cellgen(@(i) length(pts{i}), nP );

tkern=cell(2,1);
%On définit un noyau de découpe de loi

tkern{2}=@(x) (0.1+x.^2).*exp(-x.^2);
tkern{1}=@(x) exp(-x.^2);
%kernel=@(v,x) (sin(2*pi*(x-v(1)))./(x-v(1))).^2;
%tkern{2}=@(x) (x).^2 .* exp(-(x.^2));

%On choisit des moments cibles
ATmoments{1}={[-0.5 1 ;0.5 1], [1.5 4.5;2.5 7.5], [1.5 4.5;2.5 11.5]};
ATmoments{2}={[-0.5 0.5 ;0.5 1.5], [1.5 2.75;2.5 9.25],[1.5 8;2.5 8] };

AKstart{1}={[-0.5 1 ;0.5 1], [1.7 1;2.3 1],  [1.7 1;2.3 1]};
AKstart{2}={[-0.5 1 ;0.5 1], [1.7 1;2.3 1], [1.7 1;2.3 1]};

for kern=1:2
for tmom=1:2
kernel=@(v,x) tkern{kern}((x-v(1))./v(2));

% On choisit des moments cibles
Tmoments= ATmoments{tmom};






% Tmoments{k}(l,q) doit correspondre au q-ème moment de la sous-loi l de la loi parente k.
% Attention! La somme des moments partiels doit rester égale au moment de la loi parente.
% De plus, les bornes présentés, ci-haut sont très proches des bornes maximales. Les dépasser empêchera
% probablement Shred...Wishfully de converger. 


%En première approximation, on considère que les centres des noyaux devraient être proche des moments choisi
Kstart=AKstart{tmom};


% À partir de ces données ShredWithKernel construit les nL{i} sous lois
% associés aux nP lois parentes.
% Laws{i}{k} correspond à la k-ème sous-lois associés à la i-ème loi
% parente

Laws=cellgen( @(i) ShredWithKernel(Parents{i} , Tmoments{i}, Ws{i}, pts{i}, kernel, Kstart{i} ), nP) ;

for i=1:nP
path= sprintf('figs/V3/updf%d_%d_kern%d.eps',i,tmom,kern);    
divisionPlot(Parents{i},Laws{i},Ws{i},pplot,printf(path),figcount); % figcount=figcount+1;
end

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
    Mc{i}{k, mod(k,d)+1}=Laws{i}{k};
end
end

MP= @(i) Mc{i} ;

ps=[0.1,0.75];

for p=ps

q=1-p;

% Matrice de structure E
E=p*Id+q*Jd;


%% Définition de la loi matricielles
Law=matrixLawNS(A,E,MP,nP);

%% Paramètre de plot
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 


%% pdf multivarié theorique

nx=50;
ny=50;
 
%stepx= (xmax - xmin) /(nx-1);
%stepy= (ymax - ymin) /(ny-1);
 

%vy=ymin:stepy:ymax;
%vx=xmin:stepx:xmax;

vy=0:0.1:6;
vx=-3:0.1:3;

nx=length(vx);
ny=length(vy);


 [Gx,Gy]=meshgrid(vx,vy);
 
 Z=zeros(ny,nx);
 
name=@(a,b) sprintf('figs/V3/mpdf_%d-%d_%1.1f_%d_kern%d.eps',a,b,p,tmom,kern);
printingBv= @(a,b) printing(name(a,b));
 
 for i=1:nx
     for j=1:ny
        Z(j,i)=(Law.mvPartialPdf([1 2],[vx(i) vy(j)])) ;
     end
 end
figure(figcount); clf;
 contourf(Gx,Gy,Z,20);
 printingBv(1,2);
% figcount=figcount+1;
 
  for i=1:nx
     for j=1:ny
        Z(j,i)=(Law.mvPartialPdf([1 3],[vx(i) vy(j)])) ;
     end
  end
figure(figcount); clf;
 contourf(Gx,Gy,Z,20);
  printingBv(1,3);
% figcount=figcount+1;
 
  [Gx,Gy]=meshgrid(vy,vy);
  Z=zeros(ny,ny);
   for i=1:ny
     for j=1:ny
        Z(j,i)=(Law.mvPartialPdf([2 3],[vy(i) vy(j)])) ;
     end
  end
figure(figcount); clf;
 contourf(Gx,Gy,Z,20);
   printingBv(2,3);
% figcount=figcount+1;
 
 
  %% Corr
  covTheor=Law.mvMoments([1 2], [1 1]) - Law.moments(1,1)*Law.moments(1,2)
  corrTheo=covTheor/sqrt( (Law.moments(2,1)-Law.moments(1,1)^2)*(Law.moments(2,2)-Law.moments(1,2)^2) )
 
  
  %% Square Correlation
  cov2Theor=Law.mvMoments([1 2], [2 2]) - Law.moments(2,1)*Law.moments(2,2)
  corr2Theo=cov2Theor/sqrt( (Law.moments(4,1)-Law.moments(2,1)^2)*(Law.moments(4,2)-Law.moments(2,2)^2) )
 
end
end
end