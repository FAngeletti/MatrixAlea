function y = ProjHyperPlan(x,alphasum) ; 
    y = x ; 
    stmp = sum(x) ; 
    d = length(x) ; 
    if stmp ~= alphasum
        y = x + (1-stmp)/d*ones(size(x)) ; 
    end
end