addpath /Users/patriceabry/MATLAB/UTILS_STAT/

% Correlation nulle
% Correlation d'ordre 2 distincte.

clear all
close all

% taille de la série temporelle
n=10000 ;

%Matrice de structure
p=0.98;
q=1-p ;

%Loi de base
L=lnormal(0,1);

%Points de contrôles
pts=-5:0.01:5;

%Poids
Ws=[1;1]

%%Fonctions d'affinités
%f=@(x) exp(-0.25*(x.*x));
%g=@(x) exp( -(x.^2-2.5).^2);
%funAff={f g};
%% Définition des lois sectionnées
%Ls= ShredByAffinity(L, pts, Ws, funAff);

%Noyau
kernel=@(v,x) exp(-abs((x-v(1))./v(2)));
%Moments ciblés
tmoments=[0 0.2; 0 1.8 ];
% Point de départ des paramètres du noyau
kstart=tmoments;
% Définition des lois sectionnées
Ls=ShredWithKernelWishfully(L , tmoments, Ws , pts, kernel, kstart);




Lm=Ls{1};
Lp=Ls{2};

%Dimension de la matrice
d=6;

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
for(i=1:d)
Id(i,i)=1;
end 

%matrice de structure 
E=p*Id+q*J;

L2=cell(d,d);
L1=cell(d,d);


diag1={Lp,Lm,Lp,Lm,Lp,Lm};
diag2={Lp,Lp,Lp,Lm,Lm,Lm};

%Matrice de lois
L2=cell(d,d);

for i=1:d
		L2{i,i}=diag2{i};
     		L2{i,1+mod(i,d)}=diag2{i};
		L1{i,i}=diag1{i};
		L1{i,1+mod(i,d)}=diag1{i};

end

Law{2}=matrixLaw(A,E,L2,n);
Law{1}=matrixLaw(A,E,L1,n);

%lambdas=CircularEigenValues(E);
%taus=CircularTaus(lambdas)
%coeffs1= CircularEigenCoeffs(Law1.matMq(2))
%coeffs2= CircularEigenCoeffs(Law2.matMq(2))

nfig=1;
for(k=1:2)
x=Law{k}.rv();

figure(nfig); clf; nfig=nfig+1;
plot(x,'k');grid on  
set(gca,'FontSize',fontsize) 
h(1) = ylabel('Data') ; 
h(2) = title('X') ; 
set(h,'FontSize',fontsize) ; 

oplot= @(x,y,s) plot(x,y,s,'LineWidth',linewidth,'MarkerSize',markersize);
splot= @(x,s) plot(x,s,'LineWidth',linewidth,'MarkerSize',markersize);

%%  pdf

[hh,bh,gh1]=hist1d(x,50);
pdftheo=Law{k}.pdf(bh);
% figure(3) ; clf 
%   plot(bh,hh,'k'); hold on ; grid on ; 
%   plot(bh,pdftheo,'--') ; 
  
  figure(nfig) ; clf; nfig=nfig+1; 
  oplot(bh,pdftheo(:,1),'r');grid on; hold on;
  oplot(bh,hh,'k') ;
  oplot(bh, Lm.pdf(bh)./2, 'b');
  oplot(bh, Lp.pdf(bh)./2, 'b');
  
set(gca,'FontSize',fontsize) 
axis([-7 7 0 0.55])
h(1) = ylabel('Marginal') ; 
set(h,'FontSize',fontsize) ; 

% plot(bh,pdftheo(:,1),'r--','LineWidth',linewidth,'MarkerSize',markersize);
% axis(V) 
% h(1) = xlabel('lag \tau') ; 
% h(2) = ylabel('Cov X^2(\tau) X^2(0)') ; % h(3) = title(name) ; 
% set(h,'FontSize',fontsize) ; 
% set(gca,'FontSize',fontsize) 
  


 
  %% Corr
xc = xcov(x,'coeff') ;
xc=xc((n):end);


ytheo=arraygen(@(i) Law{k}.corr(i), [len]);


 
% figure(4) ; clf 
% plot(xc(1:len) ) ;  grid on; hold on;
% plot(ytheo1,'r--');
% plot(xc2(1:len),'k' ) ;  
% plot(ytheo2,'g--');
% axis(V) 

figure(nfig); clf; nfig=nfig+1;
splot(xc(1:len),'k-') ;  grid on; hold on;
splot(ytheo,'b--');
axis(V) 
%h(1) = xlabel('lag \tau') ; 
h(1) = ylabel('Corr X(\tau) X(0)') ; % h(3) = title(name) ; 
set(h,'FontSize',fontsize) ; 
set(gca,'FontSize',fontsize) ;


xcq = xcov(x.^2,'coeff') ;
xcq=xcq((n):end);



yqtheo=zeros(len,1);
v=Law{k}.genMoments([4],[1])-Law{k}.genMoments([2],[1])^2;


mq=Law{k}.moments(2);


yqtheo(1)=Law{k}.genMoments([4],[1]);



for i=2:len
yqtheo(i)=(Law{k}.genMoments([2 2],[1 (i)])-mq^2)/v;
end
 
figure(nfig); clf; nfig=nfig+1; 
splot(xcq(1:len),'k-' ) ;  grid on; hold on;
splot(yqtheo,'b--');
axis(V) 
h(1) = xlabel('lag \tau') ; 
h(2) = ylabel('Corr X^2(\tau) X^2(0)') ; % h(3) = title(name) ; 
set(h,'FontSize',fontsize) ; 
set(gca,'FontSize',fontsize) 
end


% stato(imn,channel).LWT.j1:logstato(imn,channel).LWT.j2),'-','LineWidth',linewidth,'MarkerSize',markersize,'Color',[ p p p]) ;
% set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1]) ; 
% set(gca,'XTick',[1:1:jj]) ; 
% set(gca,'XTickLabel',2.^[1:1:jj]) ;

printfig = 1 ; 
%printpath = '/Users/patriceabry/TEXTES/PUBLIS/COURANT/12PUBLIS/12ICASSP/12PRODUCT/FIGTMP/' ;
 printpath ='figs/' ; 
name = 'ExB' ; 

suffixes={ '_TS' ; '_PDF'; '_Cov';'_SqCov' };
l=length(suffixes);
if printfig == 1
for i=1:l 
    fname= @(k) sprintf('%s%s%s%d.eps',printpath,name,suffixes{i},k);
    figure(i) ; print ('-depsc', fname(1)) ;
    figure(i+l); print('-depsc', fname(2));
end
end

%  %% 
% xc = xcov(x1.^2,'coeff') ;
% xc=xc((n+1):end);
% 
% xc2 = xcov(x2.^2,'coeff') ;
% xc2=xc2((n+1):end);
% 
% len=300;
% ytheo1=zeros(len,1);
% ytheo2=zeros(len,1);
% 
% v1=Law1.genMoments([4],[1])-Law1.genMoments([2],[1])^2;
% v2=Law2.genMoments([4],[1])-Law2.genMoments([2],[1])^2;
% 
% mq1=Law1.moments(2);
% mq2=Law2.moments(2);
% 
% for i=1:len
% ytheo1(i)=(Law1.genMoments([2 2],[1 (i+1)])-mq1^2)/v1;
% ytheo2(i)=(Law2.genMoments([2 2],[1 (i+1)])-mq2^2)/v2;
% end
%  
% figure(4) ; clf 
% plot(xc(1:len) ) ;  grid on; hold on;
% plot(ytheo1,'r--');
% plot(xc2(1:len),'k' ) ;  
% plot(ytheo2,'g--');
