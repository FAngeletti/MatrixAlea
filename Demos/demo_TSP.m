
%% Version 4 + 
%% Automatic construction of the alpha parameters
% from the time scales


% size of the time series
n = 40000 ;
% size of the displayed size series
nr= 10000;








%Marginal lax
Lmarg=lLaplace(1);

%Matrix dimension of E
d=6;


%Projection operator (L(M) = <A,M>) 
A=ones(d);

% Helper matrix
J=zeros(d);
for(i=1:(d-1))
J(i,i+1 )=1;
end 
J(d,1)=1;

Id=zeros(d);
for(i=1:d)
Id(i,i)=1;
end 




%Proposed correlation lenghts
taus=[200 100 50];

lambdas= exp(-1./taus)
if (mod(d,2) ==0)
lambdas = [1 lambdas conj(lambdas((end-1):-1:1))];
else
lambdas = [1 lambdas conj(flipr(lambdas))];
end

% Alphas parameters
alphas=AlphaFromLambda(lambdas);



lambdasr=fft(alphas)

% Constructed correlation  lenghts
tausr= 1./log(abs(fft(alphas)))


% Structure matrix (polynomial decomposition)
E= alphas(d)*Id;
for i=1:(d-1)
E= J*E+ alphas(d-i)*Id;
end
E







% Construction of the matrix pdf P
%Control points
pts=-10:0.04:10;
%Weights
Ws=ones(d,1);

%Kernel
% kernel=@(v,x) exp(-abs((x-v(1))./v(2)));
kernel=@(v,x) exp(-(x-v(1)).^2./v(2));






%% Tuning of the correlation 
% Average of the marginal distribution
mu=Lmarg.moments(1);

%Correlation amplitude
c1=0.1;

%Matrix Dq(1)
D{1} = arraygen(@(n) sqrt(2*c1)*cos(pi*n), [d]);

% Tuning of the square
% Order 2 moment of the marginal distribution 
mom2=Lmarg.moments(2)

%Square correlation amplitude
c2=1.5*mom2
c3=3.5*mom2
%Matrice Dq(2) for X
D{2}{1} = arraygen(@(n) sqrt(c2)*cos(2*pi*n/d), [d])
%Matrice Dq(2) for Y
D{2}{2} = arraygen(@(n) sqrt(c3)*cos(4*pi*n/d), [d])


for k=1:2


	tmoments{k}= [(D{1}+mu)' (D{2}{k}+mom2)']

	% Initial kernel parameters
	kstart=[zeros(d,1) ones(d,1)];
         %kstart=sqrt(tmoments{k})
	% Decomposition of the marginal law in sub-laws
	Ls{k}=ShredWithKernel(Lmarg , tmoments{k}, Ws , pts, kernel, kstart);

end



L{2}=cell(d,d);
L{1}=cell(d,d);

%Law matrix
for k=1:2
for i=1:d
for diag=0:(d-1)
		L{k}{i,1+mod(i-1+diag,d)}=Ls{k}{i};
end
end
Law{k}= matrixLaw(1,E,L{k},n);
end



% display
len = 300 ; 
V = [0 len -0.1 0.5] ; 
H = [-6 6 0 0.5] ; 

fontsize = 18 ;  fontsize2 = 15 ; 
markersize = 5 ; markersize2 = 5 ; 
linewidth = 2 ; 

nfig=1;
colors={'r','b'};



%Printing options
 printpath ='figs/' ; 
name = 'ExB' ; 



suffixes={ '_TS' ; '_PDF'; '_Cov';'_SqCov' };

fname= @(suff, k) sprintf('%s%s%s%d.eps',printpath,name,suff,k);
printfig =@(i,k)   print ('-depsc', fname(suffixes{i},k))  ;


for(k=1:2)
x=Law{k}.rv();

clf;
plot(x(1:nr),'k');grid on  
set(gca,'FontSize',fontsize) 
%h(1) = ylabel('Data') ; 
%h(2) = title('X') ; 
set(gca,'FontSize',fontsize) ; 
printfig(1,k);

oplot= @(x,y,s) plot(x,y,s,'LineWidth',linewidth,'MarkerSize',markersize);
splot= @(x,s) plot(x,s,'LineWidth',linewidth,'MarkerSize',markersize);

%%  pdf

[hh,bh,gh1]=hist1d(x,100);
pdftheo=Law{k}.pdf(bh);
  solid=colors{k};
  dash=[solid '--'];

  clf;  set(gca,'FontSize',fontsize) 
  oplot(bh,pdftheo(:,1), 'k--');  axis(H); grid on;   printfig(2,k+2);

  clf;  set(gca,'FontSize',fontsize) 
  oplot(bh,pdftheo(:,1), dash);grid on; hold on;
  oplot(bh,hh,'k') ;
  for i=1:d
  oplot(bh, Ls{k}{i}.pdf(bh)./d, solid);
  end
  axis(H);
  printfig(2,k);
  

  

 
  %% Corr
xc = xcov(x,'coeff') ;
xc=xc((n):end);


ytheo=arraygen(@(i) Law{k}.corr(i), [len]);
clf;  set(gca,'FontSize',fontsize) ; 
splot(ytheo,'k--'); printfig(3,k+2); axis(V) ; grid on;
printfig(3,k+2);

clf;  set(gca,'FontSize',fontsize) ; 
splot(ytheo,dash); 
hold on;
splot(xc(1:len),'k-') ;  axis(V) ; grid on;
printfig(3,k);  




xcq = xcov(x.^2,'coeff') ;
xcq=xcq((n):end);



yqtheo=zeros(len,1);
v=Law{k}.genMoments([4],[1])-Law{k}.genMoments([2],[1])^2;


mq=Law{k}.moments(2);


yqtheo(1)=Law{k}.genMoments([4],[1]);



for i=2:len
yqtheo(i)=(Law{k}.genMoments([2 2],[1 (i)])-mq^2)/v;
end
 
clf;  set(gca,'FontSize',fontsize) ; 
splot(yqtheo,'k--'); axis(V); grid on;  printfig(4,k+2);

clf;  set(gca,'FontSize',fontsize) ; 
splot(yqtheo,dash);
hold on;
splot(xcq(1:len),'k-' ) ;   
axis(V); grid on;
printfig(4,k);
end

