# Algorithm notes

## IRLB core

The real truncated SVD kernel is a direct modern-Fortran translation of the
upstream C fast path in `src/irlb.c`. It constructs an upper bidiagonal Lanczos
matrix B, computes its small SVD with LAPACK, estimates residuals from the last
row of the left singular-vector matrix of B, and performs a thick restart with
converged Ritz vectors.

A translation trap caught during validation was the flat C indexing of the
bidiagonal superdiagonal. Because the C arrays are column-major, the expression
corresponding to the `R` coefficient is B(j,j+1), not B(j+1,j). The permanent
full-SVD regression tests specifically exercise restart iterations so that this
orientation cannot silently regress.

## Implicit transforms

For a right vector v, the effective real operator is

    (A * diag(1/scale) - 1 * center^T * diag(1/scale)) v

when both `center` and `scale` are supplied. A scalar `shift` for square
matrices is included without explicitly forming a shifted matrix. The
transpose operation applies the algebraically matching transformations.

## Sparse matrices

`csc_matrix` uses 1-based compressed-column arrays suitable for Fortran. IRLB
uses sparse matrix-vector products directly. `svdr` also uses sparse block
products directly and does not densify the input.

## Other algorithms

- `partial_eigen` follows upstream: obtain dominant singular triplets, detect a
  negative symmetric eigen-direction, shift the matrix, rerun IRLB, and remove
  the shift.
- `prcomp_irlba` performs PCA through implicitly centered/scaled IRLB.
- `ssvd` follows the Shen-Huang soft-threshold / QR iteration used upstream.
- `svdr` follows the Halko-Martinsson-Tropp block randomized iteration in
  upstream `R/svdr.R`.

## Numerical library use

Only small dense factorizations and documented fallback paths use LAPACK.
`DGESDD`, `DGEQRF`, `DORGQR`, `DSYEV`, and `ZGESDD` all have explicit Fortran
interfaces. FPM therefore does not need `implicit-external = true`.
