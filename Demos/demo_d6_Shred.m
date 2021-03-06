

clear all
close all

% Size of the time series
n=10000 ;

%Structure matrix parameters
p=0.98;
q=1-p ;

%Marginal law
 L=lLaplace(1/sqrt(2));

%Control points
pts=-5:0.01:5;

% Weights
Ws=[1;1]


%Kernel
kernel=@(v,x) exp(-abs((x-v(1))./v(2)));

%Targeted moments
tmoments=[0 0.2; 0 1.8 ];

% Initial parameters of the kernel
kstart=tmoments

% Decomposition of the marginal law in sub-laws
Ls=ShredWithKernel(L , tmoments, Ws , pts, kernel, kstart);


Lm=Ls{1};
Lp=Ls{2};

%Matrix dimension
d=6;

% Display
len = 300 ; 
V = [0 len -0.1 1] ; 
fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 

%Projection operator ( L(M) = <A,M> ) 
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

%Structure matrix
E=p*Id+q*J;

L2=cell(d,d);
L1=cell(d,d);

% D(1) : Diagonal of the moment matrix M(1)
diag1={Lp,Lm,Lp,Lm,Lp,Lm};
diag2={Lp,Lp,Lp,Lm,Lm,Lm};

%Law matrix
L2=cell(d,d);

for i=1:d
		L2{i,i}=diag2{i};
     		L2{i,1+mod(i,d)}=diag2{i};
		L1{i,i}=diag1{i};
		L1{i,1+mod(i,d)}=diag1{i};

end

% Generation of the law with matrix representation
Law2=matrixLaw(A,E,L2,n);
Law1=matrixLaw(A,E,L1,n);


x1=Law1.rv();
x2=Law2.rv();

figure(1); clf;
plot(x1,'k');grid on  
set(gca,'FontSize',fontsize) 
h(1) = ylabel('Data') ; 
h(2) = title('X') ; 
set(h,'FontSize',fontsize) ; 

figure(2); clf;
plot(x2,'k'); grid on 
set(gca,'FontSize',fontsize) 
% h(1) = ylabel('Y') ; 
h(1) = title('Y') ; 
set(h,'FontSize',fontsize) ; 

%%  pdf

[hh,bh,gh1]=hist1d(x1,50);
pdftheo=Law1.pdf(bh);

  
  figure(3) ; clf 
  plot(bh,pdftheo(:,1),'b--','LineWidth',linewidth,'MarkerSize',markersize);grid on; hold on;
plot(bh,hh,'k','LineWidth',linewidth,'MarkerSize',markersize ) ;  
set(gca,'FontSize',fontsize) 
axis([-7 7 0 0.55])
h(1) = ylabel('Marginal') ; 
set(h,'FontSize',fontsize) ; 


  
[hh,bh,gh1]=hist1d(x2,100);
pdftheo=Law2.pdf(bh);


figure(4) ; clf 
plot(bh,pdftheo(:,1),'r--','LineWidth',linewidth,'MarkerSize',markersize);grid on; hold on;
plot(bh,hh,'k','LineWidth',linewidth,'MarkerSize',markersize ) ;  
set(gca,'FontSize',fontsize) 
axis([-7 7 0 0.55])


 
  %% Corr
xc = xcov(x1,'coeff') ;
xc=xc((n):end);

xc2 = xcov(x2,'coeff') ;
xc2=xc2((n):end);

ytheo1=zeros(len,1);
ytheo2=zeros(len,1);

for i=1:len
ytheo1(i)=Law1.corr(i);
ytheo2(i)=Law2.corr(i);
end
 


figure(5) ; clf 
plot(xc(1:len),'k-','LineWidth',linewidth,'MarkerSize',markersize ) ;  grid on; hold on;
plot(ytheo1,'b--','LineWidth',linewidth,'MarkerSize',markersize);
axis(V) 
%h(1) = xlabel('lag \tau') ; 
h(1) = ylabel('Corr X(\tau) X(0)') ; % h(3) = title(name) ; 
set(h,'FontSize',fontsize) ; 
set(gca,'FontSize',fontsize) ;


figure(6); clf
plot(xc2(1:len),'k-' ,'LineWidth',linewidth,'MarkerSize',markersize) ;  grid on; hold on;
plot(ytheo2,'r--','LineWidth',linewidth,'MarkerSize',markersize);
axis(V) 
% h(1) = xlabel('lag \tau') ; 
% h(2) = ylabel('Corr X(\tau) X(0)') ; % h(3) = title(name) ; 
% set(h,'FontSize',fontsize) ; 
set(gca,'FontSize',fontsize) 

xcq = xcov(x1.^2,'coeff') ;
xcq=xcq((n):end);

xcq2 = xcov(x2.^2,'coeff') ;
xcq2=xcq2((n):end);

yqtheo1=zeros(len,1);
yqtheo2=zeros(len,1);

v1=Law1.genMoments([4],[1])-Law1.genMoments([2],[1])^2;
v2=Law2.genMoments([4],[1])-Law2.genMoments([2],[1])^2;

mq1=Law1.moments(2);
mq2=Law2.moments(2);

yqtheo1(1)=Law1.genMoments([4],[1]);
yqtheo2(1)=Law2.genMoments([4],[1]);


for i=2:len
yqtheo1(i)=(Law1.genMoments([2 2],[1 (i)])-mq1^2)/v1;
yqtheo2(i)=(Law2.genMoments([2 2],[1 (i)])-mq2^2)/v2;
end
 
figure(7) ; clf 
plot(xcq(1:len),'k-','LineWidth',linewidth,'MarkerSize',markersize ) ;  grid on; hold on;
plot(yqtheo1,'b--','LineWidth',linewidth,'MarkerSize',markersize);
axis(V) 
h(1) = xlabel('lag \tau') ; 
h(2) = ylabel('Corr X^2(\tau) X^2(0)') ; % h(3) = title(name) ; 
set(h,'FontSize',fontsize) ; 
set(gca,'FontSize',fontsize) 

figure(8); clf
plot(xcq2(1:len),'k-' ,'LineWidth',linewidth,'MarkerSize',markersize) ;  grid on; hold on;
plot(yqtheo2,'r--','LineWidth',linewidth,'MarkerSize',markersize);
axis(V) 
h(1) = xlabel('lag \tau') ; 
% h(2) = ylabel('Corr X^2(\tau) X^2(0)') ; % h(3) = title(name) ; 
set(h,'FontSize',fontsize) ; 
set(gca,'FontSize',fontsize) 

% stato(imn,channel).LWT.j1:logstato(imn,channel).LWT.j2),'-','LineWidth',linewidth,'MarkerSize',markersize,'Color',[ p p p]) ;

printfig = 0 ; 
 printpath ='../figs/' ; 
name = 'ExB' ; 

suffixes={ '_TS1' ; '_TS2'; '_PDF1';'_PDF2'; '_Cov1';'_Cov2'; '_SqCov1'; '_SqCov2' };
suffixes=cellfun( @(s) [s '.eps'], suffixes, 'UniformOutput', false);
if printfig == 1
for i=1:length(suffixes) 
    tmp_name=[printpath name suffixes{i}];
    figure(i) ; print ('-depsc', tmp_name) ; 
end
end

