function coeffs=CircularEigenCoeffs(M)
M;
L=sum(M,2);
C=sum(M,1);
n=length(L);
nvp=n-1;



% Calcul des coefficients associés à chaque valeurs propres
Lt=fft(L)./n
Cts=(fft(C))./n

coeffs=conj(Lt' .* Cts);
coeffs=coeffs(2:end);

end


