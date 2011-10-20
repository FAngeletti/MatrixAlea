function cov=CircularEigenCov(lambdas,coeffs, t)
d=length(coeffs)+1;


tmp= coeffs  .* (lambdas.^(t-1));
cov=sum(tmp);  
%if (abs(imag(cov)) >10*eps ) 
%	warning 'CircularCov : Unexpected non-zero imaginary part'
%end
cov=real(cov);
end
