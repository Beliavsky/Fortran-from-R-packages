# Lmoments-fortran

Modern Fortran/FPM translation of the computational code in the R package **Lmoments 1.3-2**.

## Scope

The library translates all substantive exports from the package NAMESPACE:

- ordinary sample L-moments and L-moment coefficients/ratios
- covariance matrices of sample L-moments
- T1 (trimmed 1,1) L-moments
- shifted Legendre coefficient matrices
- four- and six-parameter normal-polynomial quantile mixtures
- four-parameter Cauchy-polynomial quantile mixtures
- density, CDF, quantile and RNG functions for both quantile-mixture families
- parameter conversion from L-moments/T1-L-moments
- covariance transformation for the four-parameter normal-polynomial model

The API is array-based rather than R-object based. R-specific NA dispatch, matrix dimnames, and object/list formatting are intentionally not reproduced.

## Use

```fortran
use lmoments, only: dp, lmoments_sample, lmom2normpoly4, qnormpoly

real(dp) :: x(100), lm(4), param(4)
integer :: info

call lmoments_sample(x, lm, info)
param = lmom2normpoly4(lm)
print *, qnormpoly(0.95_dp, param)
```

Build/test with FPM:

```text
fpm test
fpm run --example basic
```

## Original Hosking Fortran

The user-supplied `lmoments.f` (LMOMENTS 3.04, J. R. M. Hosking, July 2005) is preserved under `upstream/`. The modernized `hosking_lmoments` reference path follows its unbiased probability-weighted-moment algorithm, while the main `lmoments_sample` routine follows the current Rcpp implementation. Hosking's AS241 normal quantile algorithm is used for the normal-polynomial distribution.

## Licensing

The R package is GPL-2.0-only. The attached Hosking library has a separate permissive IBM permission notice. The combined Fortran project is distributed under GPL-2.0-only; Hosking-derived portions retain the IBM attribution/permission notice. See `LICENSES.md` and `licenses/`.
