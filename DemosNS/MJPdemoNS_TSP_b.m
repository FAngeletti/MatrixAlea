
figcount=1;
%% Définition des lois sources
nP=2;
%Parents={lgamma(1,3), lgamma(1,3)} ;
Parents={lnormal(0,1), lgamma(2,1) };
figcount=1;

%% Découpage des sous lois
%Nombre de morceaux lois
nL={2 2};


% Ws{i}(k) Poids associé à la k-ème sous-lois de la loi source Parents{i}
Ws=cellgen(@(i) ones(nL{i},1),nP );



%Choix des points de contrôles
pts={(-4:0.01:4) (0:0.01:10) };
%On récupère le nombre de points de contrôle
npts=cellgen(@(i) length(pts{i}), 2 );


kernels{1}=@(v,x) exp(-(x-v(1)).^2./(2*v(2)));
kernels{2}=@(v,x) 1./((x-v(1)).^2+v(2)) ;
%kernel=@(v,x) (sin(2*pi*(x-v(1)))./(x-v(1))).^2;

%On choisit des moments cibles
ATmoments{1}={[-0.5 1; 
                0.5 1],  [1.5 4.5;
                          2.5 7.5]};
                      
ATmoments{2}={[-0.5 0.5 ;
                0.5 1.5], [1.5 2.75;
                           2.5 9.25]};

AKstart{1}={[-0.5 1 ;0.5 1], [1.7 1;2.3 1]};
AKstart{2}={[-0.5 0.8 ;0.5 1.2], [1.5 1;2.5 1]};
figcount=1;
for kern=1:2
    kernel=kernels{kern};
for ik=1:2
Tmoments=ATmoments{ik};    


%
% Tmoments{k}(l,q) doit correspondre au q-ème moment de la sous-loi l de la loi parente k.
% Attention! La somme des moments partiels doit rester égale au moment de la loi parente.
% De plus, les bornes présentés, ci-haut sont très proches des bornes maximales. Les dépasser empêchera
% probablement Shred...Wishfully de converger. 


%En première approximation, on considère que les centres des noyaux devraient être proche des moments choisi
Kstart=AKstart{ik};


% À partir de ces données ShredWithKernelWishfully construit les nL{i} sous lois
% associés aux nP lois parentes.
% Laws{i}{k} correspond à la k-ème sous-lois associés à la i-ème loi
% parente

Laws=cellgen( @(i) ShredWithKernelWishfully(Parents{i} , Tmoments{i}, Ws{i}, pts{i}, kernel, Kstart{i} ), nP) ;

printf=@(s) @() print('-depsc', s);
printing= @(s) print('-depsc', s);
for i=1:nP
path= sprintf('figs/V2/updf%d_%d_kern%d.eps',i,ik,kern);    
divisionPlot(Parents{i},Laws{i},Ws{i},@plot,printf(path),figcount); figcount=figcount+1;
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

ps=[0.9,0.1];

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
 

 
 for i=1:nx
     for j=1:ny
        Z(j,i)=(Law.mvPdf([vx(i) vy(j)])) ;
     end
 end
figure(figcount); clf;
 contourf(Gx,Gy,Z,20);
 name=sprintf('figs/V2/mpdf_%1.1f_%d_kern%d.eps',p,ik,kern);
 printing(name);
 figcount=figcount+1;
 
  %% Corr
  covTheor=Law.mvMoments([1 2], [1 1]) - Law.moments(1,1)*Law.moments(1,2)
  corrTheo=covTheor/sqrt( (Law.moments(2,1)-Law.moments(1,1)^2)*(Law.moments(2,2)-Law.moments(1,2)^2) )
 
  
  %% Square Correlation
  cov2Theor=Law.mvMoments([1 2], [2 2]) - Law.moments(2,1)*Law.moments(2,2)
  corr2Theo=cov2Theor/sqrt( (Law.moments(4,1)-Law.moments(2,1)^2)*(Law.moments(4,2)-Law.moments(2,2)^2) )
 
end
end
end