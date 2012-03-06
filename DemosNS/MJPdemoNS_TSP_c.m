
%% Définition des lois sources
nP=2;
%Parents={lgamma(1,3), lgamma(1,3)} ;
Parents={lnormal(0,1), lgamma(2,1) };

cellToFun=@(c) @(k) c{k};
fUniv=cellToFun(Parents);
figcount=1;

%Choix des points de contrôles
pts={(-4:0.01:4) (0:0.01:10) };
fpts=cellToFun(pts);

%On récupère le nombre de points de contrôle
npts=cellgen(@(i) length(pts{i}), 2 );


kernel=@(v,x) exp(-(x-v(1)).^2./(2*v(2)));
%kernel=@(v,x) (sin(2*pi*(x-v(1)))./(x-v(1))).^2;





%% Définition des paramètres sur la loi matricielle 1
% Moments d'ordre 1 au temps 1
M(1,:,:)= [ -0.5 -0.5; 
             0.5 0.5  ];

% Moments d'ordre 2 au temps 1
M(2,:,:) = [1 1;
            1 1];


%  Temps 2
M2(1,:,:) =[1.5 1.5; 
            2.5 2.5];
M2(2,:,:)=[4.5 4.5;
           7.5 7.5];


moments= cellToFun({M,M2});

% Paramètre de départ du noyau
pK(1,:,:)= M(1,:,:);
pK(2,:,:)= [ 1 1;
             1 1];

pK2(1,:,:)=[1.7 1.7;
            2.3 2.3] ;
pK2(2,:,:)= [ 1 1; 
              1 1];


Kstart=cellToFun({pK,pK2});


%% Création des lois matricielles
d=2;

%% Définition de A
A=ones(d);

    
%% Définition de E


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

p=0.1;
q=1-p;

% Matrice de structure E
E=p*Id+q*Jd;


%% Définition de la loi matricielles
Law=SynthesisWithConstraints(A, E , nP , fUniv, moments , kernel, Kstart, fpts ) ;
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
 figcount=figcount+1;
 
  %% Corr
  covTheor=Law.mvMoments([1 2], [1 1]) - Law.moments(1,1)*Law.moments(1,2)
  corrTheo=covTheor/sqrt( (Law.moments(2,1)-Law.moments(1,1)^2)*(Law.moments(2,2)-Law.moments(1,2)^2) )
 
  
  %% Square Correlation
  cov2Theor=Law.mvMoments([1 2], [2 2]) - Law.moments(2,1)*Law.moments(2,2)
  corr2Theo=cov2Theor/sqrt( (Law.moments(4,1)-Law.moments(2,1)^2)*(Law.moments(4,2)-Law.moments(2,2)^2) )
 
