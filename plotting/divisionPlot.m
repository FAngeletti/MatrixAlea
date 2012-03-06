function divisionPlot( Univ, Laws, Ws, plot,print, figcount )
figure(figcount); clf
Ws=Ws./sum(Ws);
F=@(y) @(x) Univ.cumulative(x)-y;

eps=0.001;
n=500;

xmin=fzero(F(eps),0); 
xmax=fzero(F(1-eps),1);

step=(xmax-xmin)/n;

v=xmin:step:xmax;

nL=length(Laws);
xlabel('x');
ylabel('p');
plot(v,Univ.pdf(v),'k'); hold on;
sty={'r--', 'b--' };
for i=1:nL
    plot(v,Ws(i).*Laws{i}.pdf(v),sty{i});
end
print();

end

