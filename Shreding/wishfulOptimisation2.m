function [ solution, err ] = wishfulOptimisation2( f, nIter,epsT, start,verbosity )
%wishfullOptimisation(f,nIter, epsT, start, verbosity) Try to optimize the function f exploiting a
%priori information on f with less than nIter steps of optimisation.
% The algorithm stops as soon as the maximum absolute error is less than epsT. 
%The verbosity varies between 0 : nothing and 4 : every steps is detailed and the variation of the maximum errer is displayed.

if(nargin<5)
	verbosity=2;
end
v=verbosity;

s=size(start);
nL=s(1);
dimM=s(2);
grad=ones(s);

if(verbosity>=4)
errv=zeros(1,nIter);
surv=zeros(1,nIter);
end

vdisp(v,1,'Wishfull optimisation algorithm starting');
solution=start;

fit=f(solution);
prev=fit;

epsD=1e-3;


for i=1:nL
for j=1:dimM
grad(i,j)=verifyGrad(f,solution,fit,epsD,i,j);
end
end

vprintf(v,'Initial gradient : \n');
vdisp(v,3,grad);

maxstep=0.5;
maxgrad= 20;

%grad= (fit-prev)./step
fit=prev;
epsR=1e-2;

sur=0.7;
nu=1.1;
relaxR=0.7;
contractR=0.5;

err=errF(fit);

vprintf(v,1,'Initial error : %f\n',err);

n=0;
while(err>epsT && n<nIter)
perr=err;
vprintf(v,2,'----- Iteration %d \n',n);
    
prev=fit;

step= -(sur*fit)./grad;
m=max(max(abs(step)));
if( m > maxstep)
step= (maxstep/m).*step;
end
vdisp(3,'Current step : '); vdisp(v,3,step); 

solution=solution+step; vdisp(v,3,'Proposed solution : '); vdisp(v,3,solution);
fit=f(solution); vdisp(v,3,'Current vectorial error : '); vdisp(v,3,fit); 

sel=find(step> epsR * m);
grad(sel)=(fit(sel)-prev(sel))./step(sel);

grad=protectGrad(f,solution,fit,grad,maxgrad,epsD);  vdisp(v,3,'Current pseudo-gradient : '); vdisp(v,3,grad); 
%grad=(fit-prev)./step;



err=errF(fit); vprintf(v,2,'Error at step %d : %f\n',n,err);

if(err<contractR*perr)
	sur=sur*nu;
elseif (perr<relaxR*err)
	sur=sur/nu;
end

n=n+1;

if(v>=4)
errv(n)=log(err);
surv(n)=log(sur);
end
vprintf(v,2,'End of iteration %d --------- \n \n',n);
end

vprintf(v,1,'Final error: %f\n',err);

if(v>=4)
figure(5);
plot(errv(1:n))

figure(6);
plot(surv(1:n))
end

end



function err=errF(fit)

err=max(max(abs(fit)));
end

function g=protectGrad(f,solution,fitn,grad,maxgrad,epsD)
g=grad;
s=size(g);
sel=find(grad<0 | grad>maxgrad );
for n=1:length(sel)
c=sel(n);
i= 1+mod(c-1,s(1));
j= 1+fix((c-1)/s(1)); 
g(c)=verifyGrad(f,solution,fitn,epsD,i,j );
end

end

function g=verifyGrad(f,solution,fitn,epsD,k,l)
tmp=solution;
tmp(k,l)= tmp(k,l)+epsD ;
tmp=f(tmp)-fitn;
g= tmp(k,l) /epsD; 
end

function vdisp(v,s,x)
	if(v>=s)
		disp(x);
	end
end

function vprintf(v,s,x,varargin)
	if(v>=s)
		fprintf(1,x,varargin{:});
	end
end

