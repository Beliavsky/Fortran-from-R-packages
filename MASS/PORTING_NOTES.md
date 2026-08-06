# Porting notes

## Representation

R model objects are replaced by derived types in `mass_types.f90` and the
vendored `rrcov_types.f90`. Optional arrays replace formula/model-frame
metadata. Status integers are returned instead of R conditions.

## Numerical equivalents

- Linear algebra is self-contained Jacobi eigenanalysis, pivoted elimination,
  and pseudoinverses rather than R's linked BLAS/LAPACK stack.
- Distribution optimization uses numerical BFGS and finite-difference Hessians.
- Robust MCD/MVE, LDA, and QDA reuse the compatible GPL-3 Fortran translation
  of `rrcov`.
- UCV/BCV/Sheather-Jones selectors evaluate direct pairwise sums instead of the
  upstream binned C implementation. The objective is the same, but floating
  point results need not be bit-for-bit identical.
- Nonmetric MDS uses pool-adjacent-violators isotonic regression and a
  deterministic gradient search.
- Stepwise AIC treats each input matrix column as one term. R formula hierarchy,
  interactions, scopes, and class-specific `extractAIC` methods are not modeled.

## Source-parity choices

The MASS `contr.sdif` coefficients and `mvrnorm(..., empirical=TRUE)` covariance
semantics are preserved. The test suite checks both explicitly.

## Licensing

MASS permits GPL version 2 or 3. The vendored robust foundation is available
under GPL-3.0-or-later. The common distributable choice for this combined work
is therefore GPL-3.0-only. Original notices and archives are preserved.
