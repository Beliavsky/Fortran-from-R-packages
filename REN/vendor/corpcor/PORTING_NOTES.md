# Porting notes

## Numerical linear algebra

The original package delegates eigendecomposition and SVD to R's numerical
libraries. This translation uses a self-contained Jacobi eigensolver for real
symmetric matrices. Compact SVD is constructed from the smaller of `A'A` and
`AA'`, retaining only singular values above the same dimension-scaled tolerance
principle used upstream.

The implementation is suitable for the covariance and shrinkage problems for
which `corpcor` is intended. Applications requiring very large dense matrices
or peak vendor-BLAS performance may replace the internal eigensolver with
LAPACK while retaining the public API.

## High-dimensional algorithms

`estimate_lambda` avoids allocating a full `p x p` matrix when `p` is much
larger than `n`. It uses the Frobenius identity
`||X'X||_F = ||XX'||_F` and computes the smaller Gram matrix.

For powers other than one, `powcor_shrink` and
`crossprod_powcor_shrink` implement the low-rank identity used by the upstream
package. Their expensive eigendecomposition is on the data rank rather than on
the full number of variables.

## R attributes

R stores shrinkage intensities and standardized partial variances as attributes.
Fortran returns them as fields of `matrix_shrinkage_result` and
`vector_shrinkage_result`.

## Missing values

The upstream routines assume a complete numeric data matrix. The Fortran port
also requires complete finite inputs. R's `NA` diagonal produced by `vec2sm`
when `diag=FALSE` is represented by IEEE quiet NaN.

## Matrix powers

Fractional powers of matrices with materially negative eigenvalues are rejected
with `corpcor_numerical_error`. Pseudoinverse powers omit numerically zero
eigenvalues, matching the shrinkage-correlation inversion behavior for
`lambda=0`.
