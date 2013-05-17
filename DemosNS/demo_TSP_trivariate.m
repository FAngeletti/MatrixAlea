
figcount=1;
printf=@(s) @() print('-depsc',s);
printing=@(s) print('-depsc',s); 


fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 3 ; 

pplot=@(x,y,s) plot(x,y,s,'LineWidth',linewidth,'MarkerSize',markersize);

%% Marginal laws
nP=3;
%Parents={lGamma(1,3), lGamma(1,3)} ;
Parents={lNormal(0,1), lGamma(2,1), lGamma(1,2) };
figcount=1;

%% Decomposition in sub-laws
%Number of sub-laws by marginal
nL={2 2 2};


% Ws{i}(k) Weights associated to the k-th sub laws of the marginal Parents{i}
Ws=cellgen(@(i) ones(nL{i},1),nP );



%Choice of the control points
pts={(-4:0.01:4) (0:0.01:10) (0:0.01:10) };
%Number of control points
npts=cellgen(@(i) length(pts{i}), nP );

% Definition of the decomposition kernels
tkern=cell(2,1);


tkern{2}=@(x) (0.1+x.^2).*exp(-x.^2);
tkern{1}=@(x) exp(-x.^2);
%kernel=@(v,x) (sin(2*pi*(x-v(1)))./(x-v(1))).^2;
%tkern{2}=@(x) (x).^2 .* exp(-(x.^2));

%Targeted moments
ATmoments{1}={[-0.5 1 ;0.5 1], [1.5 4.5;2.5 7.5], [1.5 4.5;2.5 11.5]};
ATmoments{2}={[-0.5 0.5 ;0.5 1.5], [1.5 2.75;2.5 9.25],[1.5 8;2.5 8] };

AKstart{1}={[-0.5 1 ;0.5 1], [1.7 1;2.3 1],  [1.7 1;2.3 1]};
AKstart{2}={[-0.5 1 ;0.5 1], [1.7 1;2.3 1], [1.7 1;2.3 1]};

% Construction of the law matrix
for kern=1:2
for tmom=1:2

% kernel
kernel=@(v,x) tkern{kern}((x-v(1))./v(2));

%Targeted moments
Tmoments= ATmoments{tmom};






% Tmoments{k}(l,q) must correspond to the moment of order q of the sublaw l  of the marginal  k.
% Warning! The sums of the paryial moments must stay equal to 
the moment of the marginal law. 
% Moreover, the bounds used here are probably near the maximal bounds. 
%The convergence of ShredWithKernel is very uncertain
beyond them.

 

%Intial parameters for the kernl : 
Kstart=AKstart{tmom};


% Construction of the decompositions of the marginal distributions
% Laws{i}{k} correspond to the k-th sub-laws of the i-th marginal distribution

Laws=cellgen( @(i) ShredWithKernel(Parents{i} , Tmoments{i}, Ws{i}, pts{i}, kernel, Kstart{i} ), nP) ;

for i=1:nP
% Figures are printed in a figs/V3 subdirectory
path= sprintf('figs/V3/updf%d_%d_kern%d.eps',i,tmom,kern);    
divisionPlot(Parents{i},Laws{i},Ws{i},pplot,printf(path),figcount); % figcount=figcount+1;
end

%% Dimension of E
d=2;

% Idendity matrix
Id=zeros(d);
for i=1:d
    Id(i,i)=1;
end

% Circulant matrix
Jd=zeros(d);
for i=1:d
    Jd(i, mod(i, d) +1)=1;
end


%% Définition of A
A=ones(d);

%% Construction of the law (pdf) matrix function MP(k)
% MP(k) is the matrix probability at time k 
% Mc is the cell used to stocked MP(k)
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

% Creation of MP 
MP= @(i) Mc{i} ;

% Different values of p
ps=[0.1,0.8];

for p=ps

q=1-p;

% Structure matrix E
E=p*Id+q*Jd;


%% Construction of the random vector with matrix representation
Law=matrixLawNS(A,E,MP,nP);

%% Plotting parameters
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 


%% Theoretical multivariate pdf
nx=50;
ny=50;
 


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

colormap( flipud( gray(20) ) );

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
