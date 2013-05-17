
%% Shared (non-zero) correlation 
%% Different square correlation

%Matrix dimension
d=8;


%Projection operator (L(M) = <A,M>) 
A=ones(d);



%Marginal law
L=lnormal(0,1);

%Control points
pts=-5:0.001:5;

%Weights
Ws=ones(4,1);

%Kernel
kernel=@(v,x) exp(-abs((x-v(1))./v(2)));

%Targeted high average
mh=0.5;

%Maximal high variance : vh=1.7
vh=1.6;

vl=2-vh;ml= -mh;

% Targeted moments
tmoments=[ mh vh; ml vh; mh vl; ml vl ];

%Initial kernel parameters
kstart=tmoments

% Decomposition of the marginal distribution in sub-laws
Ls=ShredWithKernel(L , tmoments, Ws , pts, kernel, kstart);

% Renaming
L0=L;
Lpsm=Ls{3};
Lpsp=Ls{1};
Lmsm=Ls{4};
Lmsp=Ls{2};



%% Shared marginal and correlation
%% Different square correlation 
% diag%d : Diagonal of the moment matrix 

diag1={Lpsp Lmsp Lpsp Lmsp Lpsm Lmsm Lpsm Lmsm};
diag2={Lpsp Lmsp Lpsm Lmsm Lpsp Lmsp Lpsm Lmsm};

%Law matrix
L1=cell(d,d);
L2=cell(d,d);


for i=1:d
    for j=1:d
	if(i==j)
		L1{i,j}=diag1{i};
		L2{i,j}=diag2{i};
	else
     		L1{i,j}=L0 ;
		L2{i,j}=L0;
	end
    end
end

%Structure matrix construction
p=0.98;
q=1-p ;

J=zeros(d);
for(i=1:(d-1))
J(i,i+1 )=1;
end 
J(d,1)=1;

Id=zeros(d);
for(i=1:d)
Id(i,i)=1;
end 

% Structure matrix
E= p.*Id+q.*(J) ;

%Size of the time series
n=500000;

Law1=matrixLaw(A, E , L1 , n);
Law2=matrixLaw(A,E,L2,n);

x1=Law1.rv();
x2=Law2.rv();

figure(1); clf;
plot(x1);

figure(2); clf;
plot(x2);






[hh,bh,gh1]=hist1d(x1,100);
pdftheo=Law1.pdf(bh);

figure(3) ; clf 
  plot(bh,hh,'k'); hold on ; grid on ; 
  plot(bh,pdftheo,'--') ; 
  
[hh,bh,gh1]=hist1d(x2,100);
pdftheo=Law2.pdf(bh);

 plot(bh,hh,'b');
 plot(bh,pdftheo,'g--') ; 


xc = xcov(x1,'coeff') ;
xc=xc((n+1):end);

xc2 = xcov(x2,'coeff') ;
xc2=xc2((n+1):end);

len=300;
ytheo1=zeros(len,1);
ytheo2=zeros(len,1);

for i=1:len
ytheo1(i)=Law1.corr(i+1);
ytheo2(i)=Law2.corr(i+1);
end
 
figure(4) ; clf 
plot(xc(1:len) ) ;  grid on; hold on;
plot(ytheo1,'r--');
plot(xc2(1:len),'k' ) ;  
plot(ytheo2,'g--');


xcq = xcov(x1.^2,'coeff') ;
xcq=xcq((n+1):end);

xcq2 = xcov(x2.^2,'coeff') ;
xcq2=xcq2((n+1):end);

len=300;
yqtheo1=zeros(len,1);
yqtheo2=zeros(len,1);

v1=Law1.genMoments([4],[1])-Law1.genMoments([2],[1])^2;
v2=Law2.genMoments([4],[1])-Law2.genMoments([2],[1])^2;

mq1=Law1.moments(2);

mq2=Law2.moments(2);

for i=1:len
yqtheo1(i)=(Law1.genMoments([2 2],[1 (i+1)])-mq1^2)/v1;
yqtheo2(i)=(Law2.genMoments([2 2],[1 (i+1)])-mq2^2)/v2;
end
 
figure(5) ; clf 
plot(xcq(1:len) ) ;  grid on; hold on;
plot(yqtheo1,'r--');
plot(xcq2(1:len),'k' ) ;  
plot(yqtheo2,'g--');
