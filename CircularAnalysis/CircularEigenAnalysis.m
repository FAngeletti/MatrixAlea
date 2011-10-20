function [lambdas,taus,Ts,coeffs]=CircularEigenAnalysis(Ml)

E=Ml.E;
M=Ml.matMq(1);

lambdas=CircularEigenValues(E);
taus=CircularTaus(lambdas);
Ts=CircularTs(lambdas);
coeffs=CircularEigenCoeffs(M);

end
