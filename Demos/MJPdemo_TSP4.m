addpath /Users/pabry/MATLAB/UTILS_STAT/
% remplace marginale gaussienne par marginale student

% Correlation nulle
% Correlation d'ordre 2 distincte.

clear all
% close all

% taille de la série temporelle
n = 40000 ;
nr= 10000;

%Matrice de structure
p=0.99;
q=1-p ;

%Loi de base
Lmarg=lLaplace(1);

%Points de contrôles
pts=-10:0.04:10;

%Dimension de la matrice E
d=6;

%Poids
Ws=ones(d,1);


%Noyau
 kernel=@(v,x) exp(-abs((x-v(1))./v(2)));
%kernel=@(v,x) exp(-(x-v(1)).^2./v(2));






% Réglage de la correlation
% Moyenne de la distribution stationnaire
mu=Lmarg.moments(1);
%Amplitude de la correlation
c1=0.1;
%Matrice Dq(1)
D{1} = arraygen(@(n) sqrt(2*c1)*cos(pi*n), [d]);

% Réglage de la correlation des carrés
% Moment d'ordre 2 de la distribution stationnaire
mom2=Lmarg.moments(2)

%Amplitude de la correlation
c2=1.5*mom2
c3=3.5*mom2
%Matrice Dq(2) pour X
D{2}{1} = arraygen(@(n) sqrt(c2)*cos(2*pi*n/d), [d])
%Matrice Dq(2) pour Y
D{2}{2} = arraygen(@(n) sqrt(c3)*cos(4*pi*n/d), [d])


for k=1:2


	tmoments{k}= [(D{1}+mu)' (D{2}{k}+mom2)']

	% Point de départ des paramètres du noyau
	kstart=tmoments{k};
	% Définition des lois sectionnées
	Ls{k}=ShredWithKernel(Lmarg , tmoments{k}, Ws , pts, kernel, kstart);

end


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

L{2}=cell(d,d);
L{1}=cell(d,d);

%Matrice des lois
for k=1:2
for i=1:d
		L{k}{i,i}=Ls{k}{i};
     		L{k}{i,1+mod(i,d)}=Ls{k}{i};
end
Law{k}= matrixLaw(1,E,L{k},n);
end


%lambdas=CircularEigenValues(E);
%taus=CircularTaus(lambdas)
%coeffs1= CircularEigenCoeffs(Law1.matMq(2))
%coeffs2= CircularEigenCoeffs(Law2.matMq(2))



% affichage
len = 300 ; 
V = [0 len -0.1 1] ; 
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 

nfig=1;
colors={'r','b'};





for(k=1:2)
x=Law{k}.rv();

figure(nfig); clf; nfig=nfig+1;
plot(x(1:nr),'k');grid on  
set(gca,'FontSize',fontsize) 
h(1) = ylabel('Data') ; 
h(2) = title('X') ; 
set(h,'FontSize',fontsize) ; 

oplot= @(x,y,s) plot(x,y,s,'LineWidth',linewidth,'MarkerSize',markersize);
splot= @(x,s) plot(x,s,'LineWidth',linewidth,'MarkerSize',markersize);

%%  pdf

[hh,bh,gh1]=hist1d(x,100);
pdftheo=Law{k}.pdf(bh);
% figure(3) ; clf 
%   plot(bh,hh,'k'); hold on ; grid on ; 
%   plot(bh,pdftheo,'--') ; 
  solid=colors{k};
  dash=[solid '--'];
  figure(nfig) ; clf; nfig=nfig+1; 
  oplot(bh,pdftheo(:,1), dash);grid on; hold on;
  oplot(bh,hh,'k') ;
  for i=1:d
  oplot(bh, Ls{k}{i}.pdf(bh)./d, solid);
  end
  
set(gca,'FontSize',fontsize) 
%axis([-7 7 0 0.55])
%h(1) = ylabel('Marginal') ; 
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
splot(ytheo,dash);
axis(V) 
%h(1) = xlabel('lag \tau') ; 
%h(1) = ylabel('Corr X(\tau) X(0)') ; % h(3) = title(name) ; 
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
splot(yqtheo,dash);
axis(V) 
%h(1) = xlabel('lag \tau') ; 
%h(2) = ylabel('Corr X^2(\tau) X^2(0)') ; % h(3) = title(name) ; 
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
