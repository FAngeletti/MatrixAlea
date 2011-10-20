d=120;
c=cell(d,d);
f= @(x) (@(y) (x*y) );
for i=1:d
for j=1:d
	c{i,j}=lnormal(i,j);
end
end

x=0;
for i=1:d
for j=1:d
	l=c{i,j};
	x=x+l.rv();
end
end
x
