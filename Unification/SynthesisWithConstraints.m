function ML=SynthesisWithConstraints(A, E , n , UnivariateLaws, Moments , kernel, kernel_parameters0, pts ) 
%% Create a matrix law with linear form L= tr(A' .) , structure matrice E sastifying the required constraints
% A :matrix d-by-d
% E : matrix d-by-d
% n : size of the signal
% UnivariateLaws : n->Univariate law : function mapping k to the kth univariate law of the signal
% Moments : n-> matrix (o*d*d) where o is the number of constrained moments.
% pts : n-> array discretisation points.
% kernel : the kernel function : (P,R) -> R  where P is a parameter space
% kernel_parameter0 : n->P initial parameter of the kernel.

Ps=cell(n,1);
for k=1:n


norm = sum(sum(A.*E^n));
Ws= (E^(n-k)*A'*E^(k-1)).*E./norm ;
d=length(A);
M=Moments(k);
Ws= reshape(Ws,[ d*d 1]);
s=size(M);
m=reshape(M, [s(1) s(2)*s(3)]);
m=m';

kp=reshape(kernel_parameters0(k), [s(1) s(2)*s(3)]);
kp=kp';

Laws=ShredWithKernel(UnivariateLaws(k) , m, Ws, pts(k), kernel, kp );

p=cell(d,d);

for i=1:d
  for j=1:d
    p{i,j}=Laws{i+d*(j-1)};
  end
end
Ps{k}=p;
end
fP=@(n) Ps{n};

ML=matrixLawNS(A,E,fP,n);



end
