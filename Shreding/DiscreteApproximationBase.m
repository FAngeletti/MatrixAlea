function [ points, EpsT ] = DiscreteApproximationBase(L,eps,start,final)
%DiscreteApproximationBase(L,eps,start,end) Create an approximation base for the law L

points(1)=start;
pos=start;
step=0.1;
extr=L.extrema();
i=2;
s=0;
EpsT=1+L.cumulative(start)-L.cumulative(final);
for p=[extr final]
    
while(error(L,s,pos,p)>eps)
    f= @(x) error(L,s,pos,pos+abs(x)) -eps;
    step=abs(fzero(f,step));
    EpsT=EpsT+error(L,s,pos,pos+step);
    pos=pos+step;
    points(i)=pos; 
    i=i+1;
end
s=s-1;
end



end

function err=error(L,s,x,y)
l=abs(y-x);
if(s==0)
    p=l*L.pdf(x);
else
    p=l*L.pdf(y);
end
rP=L.cumulative(y)-L.cumulative(x);
err=rP - p;
end

