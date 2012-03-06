L=lnormal(0,1);
%L=lgamma(2,1);

pts=-4:0.01:4;
%pts=0:0.01:10;

nL=8;
Ws=ones(nL,1);
kernel=@(v,x) exp(-((x-v(1))).^2);
%kernel=@(v,y) 1./(1+((v(1)-y)).^2); 
%kernel=@(v,x) exp(-abs((x-v(1))));
%kernel=@(v,x) abs(sinc(x-v(1)));

%kernel=@(v,x) exp(-((x-v(1))./v(2)).^2);
%kernel=@(v,y) 1./(1+((v(1)-y)./v(2)).^2); 
%kernel=@(v,x) exp(-abs((x-v(1))./v(2)));



mh=1.2;
pas=2*mh/(nL-1);

tmoments=[-1; -0.75; -0.5;-0.25;0.25;0.5;0.75;1];
%tmoments=vertcat(tmomentsH,-tmomentsH);
%tmoments= (-mh:pas:mh)';
%tmoments=[-0.4 0.3;-0.4 1.35; .8  1.35]
start=tmoments;

tic;
Ls=ShredWithKernelWishfully(L , tmoments, Ws, pts, kernel,start);
toc

TestMixture(L,Ls,Ws,100,50);

for i=1:nL
    fprintf('E X^2 = %d\n' ,Ls{i}.moments(2));
end
