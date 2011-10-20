function [ y ] = InvCumulative(f, x )
%InvCumulative Inversion intelligente de cumulative

% Bound search
sup=sbound(f,x);
inf=ibound(f,x);


y=fzero(@(t) f(t)-x, [inf sup] );

end

function s=sbound(f,x)
st=1;
t=f(st);

while(t<x)
st=2*st;
t=f(st);
end

s=st;
end



function inf=ibound(f,x)
st=-1;
t=f(st);

while(t>x)
st=2*st;
t=fe(st);
end

inf=st;
end