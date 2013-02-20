function alpha=AlphaFromLambda(lambda)

d=length(lambda);
nbiterationmax = 200 ;
alphasum = 1 ; 
m = 0 ; 
M = 1 ; 
gamman = 0.49 ; % in [0, 1] 

flambda = ifft(lambda);

zn = zeros(1,d) ; 
alphanold = zn ; 
tn = 1 ; 
for n =1:1:nbiterationmax
    g = zn-2*gamman*(zn-flambda) ; 
    y = ProjRange(g,m,M) ;
    alphan = ProjHyperPlan(y,alphasum) ;
    tnnew = (1+sqrt(1+4*tn^2))/2 ; 
    zn = alphan + (tn-1)/tnnew*(alphan-alphanold) ; 
    alphanold = alphan ; 
    tn = tnnew ; 
    crit(n) = sum((alphan-flambda).^2) ; 
    dc1(n) = norm(alphan - ProjHyperPlan(alphan,alphasum)) ; 
    dc2(n) = norm(alphan - ProjRange(alphan,m,M) ) ; 
    % fprintf('\n % 3.2f\t % 3.2f\t % 3.2f\n',crit(n),dc1(n),dc2(n)) ; 
end

    alpha = ProjHyperPlan(alphan,alphasum) ;
    alpha = ProjRange(alpha,m,M) ;
    
alpha=real(alpha);
end

