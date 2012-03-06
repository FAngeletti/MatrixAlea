% Définition de la loi parente
Parent=lgamma(2,1);

%Nombre de lois
nL=2;

% Poids associé à chaque lois
Ws=ones(nL,1);

%Définitions des points de contrôle
nps=5000;
pts=vRandL(Parent,nps);
%pts=[-2,-1,0.5,0.75,1,2];

%Définition des fonctions d'affinité
kernel=@(s,x,y) 1./(1+((x-y)/s).^2); 

funAff=cell(nL,1);
for i=1:nL
   % funAff{i}=@(x) kernel(0.5,2*(2*(i-1)-(nL-1))/(nL-1),x) ;
   funAff{i}=@(x) kernel(0.5,4*(i/nL)-2,x) ;
end

Laws=ShredByAffinity(Parent,pts, Ws, funAff);

% Loi de contrôle
Remix=lMixture(Laws, Ws);

%TestLaw(Remix,1000,50, 1:1:4);
%CompareLaws({Parent,Remix}, 1000, 50, 1:1:4);

TestMixture(Parent,Laws,Ws,2000,50);

