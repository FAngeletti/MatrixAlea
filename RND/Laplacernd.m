% Laplacernd.m
% 

function x = Laplacernd(s) 

tmp = sign(rand-1/2) ; 
while tmp == 0 
    tmp = sign(rand-1/2) ; 
end
x = tmp.*exprnd(s) ; 
end