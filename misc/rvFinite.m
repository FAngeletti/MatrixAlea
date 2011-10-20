function i=rvFinite(v)
n=length(v);
p=cumsum(v);
p=p./p(n);

vmin=0;
vmax=n;

u=unifrnd(0,1);
if (u<p(1))
	i=1;
else

while(vmax-vmin>1)	
vmid=floor( (vmin+vmax)/2 );
if(p(vmid)>u)
	vmax=vmid;
else 
	vmin=vmid;
end
i=vmax;
end
end
