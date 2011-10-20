
%%%%%%%%
%Fonction de test de Loi
%
function TestLaw(L,n,nh, qs)
%%Synthesis
x=vRandL(L,n);
[h,b,g]=hist1d(x,nh);

%% PDF Test
figure(1);clf
plot(b,h,'k'); hold on; plot(b, L.pdf(b),'r--' ); 


%% Cumulative test
figure(2);clf
xs=sort(x);
fs=0.*xs;
for i=1:n
    fs(i)=L.cumulative(xs(i));
end
plot(xs, (1:n)./n,'k'); hold on; plot(xs, fs,'r--');

%% Moments test
ms=0.*qs;
tms=0.*qs;
for i=1:length(qs)
    tms(i)=(L.moments(qs(i)));
    ms(i)= (sum(x.^qs(i) )/n);
end

figure(3); clf
plot(qs,ms,'k'); hold on; plot(qs,tms,'r--')






end