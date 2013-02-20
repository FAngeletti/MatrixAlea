#Presentation 

This matlab package gathers a set of functions and subroutines for the synthesis and analysis of random variables
with matrix representation. In brief, The joint probability density of these random variables
can be expressed as 
p(x_1, ..., x_n) = L (R_n(x_1 ).. R_n(x_n) ) / Z 
where R is a positive matrix function and L is a linear form on the space of matrices.
These random variables with a matrix joint probability density enjoy various statistical
properties which are explored in a (future) series of papers (cf http://arxiv.org/abs/1203.4500).
This package aims to provide a clear API in order to helps the practitionner to explore
the potentiality of such laws. This aim is still faraway for now, so you are encouraged to provide
comments and critics.

#Structure of the package
## Laws class
In order to implement probability laws, this package use 2 main abstract class.
The univariate laws implements the abstract class lUniv found in the subdirectory LawsUniv.
A certain number of reimplementation of classical laws are also provided in this subdirectory.

Multivariate laws implements the abstract class lMulti found in the subdirectory LawsMulti.

##Random variables with matrix representation
The core of this package are the class "matrixLaw" and "matrixLawNS" found in the subdirectory LawsMulti.
 The first class corresponds to the case where R_n is constant in n, the second class models the general situation.

### Homogeneous Matrix laws
In order to create a random variable with matrix variables, one needs to use
```matlab
mL= matrixLaw(A,E,Laws,n)
```
In the previous code fragment, 
* 'A' is a positive positive matrix which corresponds to the dual vector of L, i.e. L(M)=tr( A^T E).
* 'E' is a positive matrix which corresponds to the integral of the matrix R_n(x). 
* 'Laws' is a 2-dimensionnal cell of laws.

For more details, see the tutorial/Start.m.

### Inhomogeneous matrix laws
The inhomogeneous case looks quite similar to the homogeneous case 
```matlab
mL= matrixLawNS(A,E,fLaws,n)
```
 except that fLaws must be a function such that fLaws(k) is the cells of laws at
the time k.

## Laws shreding
In typical utilisation, it is more practical to generate random variables with prescribed
marginal(s) and dependency structures rather to prescribes the whole joint probability density function.
The aims of the "Shreding" subdirectory is to provide helpers functions, which are able to shreds
a parent law into a set of sub-laws such that the parent law can be recomposed as a mixture of these sublaws.

The end user is gently but strongly encouraged to only use 'ShredWithKernel' which is the only human-friendly function available for now (non-human users should contact directly the author). 

```matlab
subLaws = ShredWithKernel(L , moments, Ws0, pts0, Kern, Kstart )
```

* 'L' is the parent Law
* 'moments' contains the constraints on the sub-laws moments, which should be a s*o matrix where s is the number of sublaws and o is the order of the moments constraint.
* 'Ws0' is a vector of dimension s, which contains the weight of the mixture, i.e. L= sum( wWs0_i * subLaws(i) )
* pts0 is a vector containing the discretization points which would be used for the constructions of the sub-laws L
* Kern is the kernel function, such that Kern(param) is a real endofunction.
* Kstart is the starting parameters for the kernel function

## Synthesis with constraints

To be written
## Synthesis with circular matrix

To be written

## Demos and example

Many (non exactly understandable) examples  are present in the Demos and DemosNS subdirectories.

