
%%%%%%%%
%Fonction de comparaison de Lois
%
function CompareLaws(Ls,n, nh, qs)

colors={'r' 'b' 'g' 'k'};
nL=length(Ls);

%%Synthesis
x=cell(1,nL);
for i=1:nL
x{i}=vRandL(Ls{i},n);
end



%% PDF Test
figure(1);clf; hold on;
figure(2);clf; hold on;

for k=1:nL
[h,b,g]=hist1d(x{k},nh);
figure(1);
plot(b,h);
figure(2);
plot(b, Ls{k}.pdf(b), colors{k} );
end

%% Cumulative test
figure(3);clf; hold on
figure(4); clf; hold on;

for k=1:nL
xs=sort(x{k});
fs=arrayfun(@(x) Ls{k}.cumulative(x), xs);

figure(3);
plot(xs, (1:n)./n,colors{k} );
figure(4); plot(xs, fs,colors{k});
end

%% Moments test
figure(4);clf; hold on
figure(5); clf; hold on;

for k=1:nL
ms=arrayfun( @(q) log(1+abs(sum(x{k}.^q )/n)), qs);
tms=arrayfun( @(q) log(1+abs(Ls{k}.moments(q))), qs);

figure(4); 
plot(qs,ms,colors{k});

figure(5)
plot(qs,tms,colors{k});

end





end