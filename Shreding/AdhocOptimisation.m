function [ solution, err ] = wishfulOptimisation( f, nIter, start )
%wishfullOptimisation(f,start) Try to optimize the function f exploiting a
%priori information on f.

s=size(start);
nL=s(1);
dimM=s(2);
grad=ones(s);


solution=start;

fit=f(solution)
prev=fit;


step= -(1e-3*fit)./grad
solution=solution+step
prev=fit;
fit=f(solution+step)

grad=(fit-prev)./step
fit=prev;

maxstep=0.25.*ones(s);


%grad= (fit-prev)./step
fit=prev;
epsT=1e-5;
epsR=1e-2;

sur=1.1;
nu=1.5;
err=errF(fit);


n=0;
while(err>epsT && n<nIter)
fprintf(1,'Iteration %d \n',n);
    
prev=fit;

step= -(sur*fit)./grad
m=max(max(step));
if( m > maxstep)
step= (maxstep/m).*step;
end
solution=solution+step
fit=f(solution)

sel=find(step> epsR * m);
grad(sel)=(fit(sel)-prev(sel))./step(sel)
%grad=(fit-prev)./step;
perr=err;
err=errF(fit)


n=n+1;

fprintf(1,'Fin de l iteration %d --------- \n \n',n);
end
n

end



function err=errF(fit)

err=max(max(abs(fit)));
end

function g=verifyGrad(solution,fit0,epsD,k,l)
solution(k,l)= solution(k,l)+epsD ;
g=fit(0); 
end

