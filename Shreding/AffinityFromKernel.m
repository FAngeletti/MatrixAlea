function funAff = AffinityFromKernel( kernel, kernelParams )
%AffinityKernel(kernel , kernelArg) Generate a cell of affinity function using a kernel function and its argument
% kernel(param, x) : kernel function with parameter argument arg at point x
% kernelParams(i, k) : nShred * dParam matrice. nShred number of shred law.
% nArg : dimension of the vector parameter of the kernel
% kernelParams(i,:) : parameter vector for affinity function at point x

funAff=cellgen( @(i) @(x) kernel( kernelParams(i,:),x), length(kernelParams(1,:)) );

end

