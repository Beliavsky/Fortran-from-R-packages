# bigstatsr-fortran 0.1.0

**Official CRAN title:** Statistical Tools for Filebacked Big Matrices

Modern Fortran computational-core translation of `bigstatsr` 1.6.2, packaged
for the Fortran Package Manager (FPM).

The upstream R package combines an R/S4 file-backed-matrix interface with
Rcpp kernels, RSpectra, parallel orchestration, plotting, and model wrappers.
This port concentrates on the numerical kernels and represents an FBM with a
portable Fortran stream-backed matrix type.

## Implemented

### File-backed matrices

- `fbm_real`: column-major real64 backing files
- `fbm_code256`: raw 0..255 backing files
- create/attach/read/write/get/set operations
- array-to-FBM and FBM-to-array conversion
- copy, increment, and blocked transpose
- code256 row/column counts

The stream-backed design does not require Rcpp or POSIX `mmap`, and therefore
works with ordinary Fortran runtimes on Windows, macOS, and Unix-like systems.
It trades some random-access performance for portability.

### Matrix/statistical kernels

- `big_colstats`: sums and sample variances
- `big_scale`: center/scale vectors
- `big_prod_vec`, `big_cprod_vec`
- `big_prod_mat`, `big_cprod_mat`
- `big_crossprod_self`, `big_tcrossprod_self`
- `big_cor`
- AUC and bootstrap AUC
- numeric-matrix partial correlation with Fisher CI
- grouped summary kernels used by sparse regression
- `block_size`, `rows_along`, `cols_along`, `get_beta`

### Regression

- column-wise linear regression with covariate adjustment
- SVD rank filtering matching `big_univLinReg`
- column-wise logistic IRLS with covariance/SEs
- Gaussian elastic-net/lasso coefficient paths
- logistic elastic-net/lasso coefficient paths
- coefficient prediction helper

The elastic-net kernels target the same objective as the C++ coordinate-descent
code in `bigstatsr`. The R package's outer CMSA cross-validation/grid-search
wrapper is not reproduced in v0.1.0.

### PCA / SVD

- primal/dual `big_SVD` analogue
- matrix-free `big_randomSVD` analogue
- centering and scaling
- SVD prediction / principal-component scores

The matrix-free path uses the vendored MPL-2.0 Fortran RSpectra translation and
the repository's shared free-form ARPACK-NG and pure-Fortran LAPACK backends.

## Example

```fortran
program demo
    use bigstatsr
    implicit none
    type(fbm_real), target :: x
    type(colstats_result) :: st
    type(big_svd_result) :: fit
    real(dp) :: a(5,3)

    a = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
                 2.0_dp,1.0_dp,0.0_dp,-1.0_dp,-2.0_dp, &
                 1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp], [5,3])

    x = create_fbm("x.bk", 5, 3)
    call fbm_from_array(x, a)
    st = big_colstats(x)
    fit = big_svd(x, 2, center=.true.)
end program demo
```

## Build

The default FPM build requires no system ARPACK, LAPACK, or BLAS library:

```text
fpm test
```

## Validation

The release is tested with GNU Fortran using:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

Tests cover file-backed I/O, matrix products, counts, AUC/partial correlation,
univariate regression, Gaussian/logistic elastic-net paths, primal/dual SVD,
and matrix-free ARPACK SVD.

## Scope

This is a computational-core translation, not an R runtime emulator. R S4 and
reference-class behavior, RDS metadata, plotting, `foreach`/cluster management,
data-frame factor encoding, `bigreadr` text import/export, and S3 presentation
methods are intentionally outside the Fortran API. See
`docs/TRANSLATION_COVERAGE.md` for the export-by-export mapping.
