function  TestMixture( Pl, Mxs, ws0 , n, nh)
%TestMixture(Pl,Mxs,ws0,n,nh) Verify that the laws Mxs mixed with weight ws recreate Pl 
% 
ws=ws0./sum(ws0);
nl=length(Mxs);
xp=vRandL(Pl,n);
xms=cell(1,nl);

Xmax=max(xp);
Xmin=min(xp);

for i=1:nl
    xms{i}=vRandL(Mxs{i},n);
    Xmax=max(Xmax,max(xms{i}));
    Xmin=min(Xmin,min(xms{i}));
end

l=Xmax-Xmin;
bh=Xmin:l/(nh+1):Xmax;
center=( bh(2:end) + bh(1:(end-1) ) )/2;

colors=hsv(nl);

figure(1);clf; hold on;
h=hist(bh,xp);
plot(center,h,'k');
hT=0.*h;
for i=1:nl
    hT=ws(i).*hist(bh,xms{i})+hT;
    plot(center,hT, 'color', colors(i,:) ); 
end

figure(2); clf; hold on;
pdfT=0.*center;

plot(center,Pl.pdf(center),'k' );

for i=1:nl
    pdfT=ws(i).*Mxs{i}.pdf(center)+pdfT;
    plot(center, pdfT,  'color', colors(i,:));
end

figure(3);clf; hold on;
plot(center,Pl.pdf(center),'k');
for i=1:nl
    pdfT=ws(i).*Mxs{i}.pdf(center);
    plot(center, pdfT, 'color', colors(i,:) );
end


end

function h=hist(bh,x)
n=length(bh)-1;
l=bh(2)-bh(1);
xorg=bh(1);
h=zeros(1,n);
np=length(x);
for y=x
    pos= fix( (y-xorg) /l);
    if(pos>0 && pos <=n )
        h(pos)=h(pos)+1;
    end
end

h=h./(np *diff(bh) ); 

end

