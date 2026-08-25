# Porting notes

## Gradient projection

The orthogonal and oblique algorithms follow the 2026.8-2 R implementation:

- projected gradient on the orthogonal/oblique constraint manifold;
- Barzilai-Borwein step length with the same 1e-10 to 20 clipping;
- non-monotone objective window;
- 11-step backtracking line search;
- column normalization for the oblique transformation;
- optional Cayley step for orthogonal rotation.

The R implementation uses `svd(X)$u %*% t(svd(X)$v)` for the orthogonal
snap-back. The Fortran implementation computes the same polar factor as
`X (X'X)^(-1/2)` using a symmetric Jacobi eigensolver. This avoids a BLAS/LAPACK
runtime dependency while representing the same nearest-orthogonal projection.

## Random starts

`random_start` uses standard-normal matrices followed by a QR-equivalent
modified Gram-Schmidt factorization with positive diagonal convention. For a
Gaussian input matrix this gives the same Haar orthogonal distribution targeted
by the sign-corrected R QR implementation.

## Oblique singularity handling

The R package falls back to an SVD pseudoinverse if a transformation becomes
singular. The standalone Fortran engine reports a nonzero `info` code instead
of silently continuing with a rank-deficient transformation. Well-conditioned
factor rotations do not use this path. This is a deliberate safety difference.

## Target matrices

R uses `NA` entries to mark unspecified target elements. In Fortran,
`vgq_target` treats IEEE NaN target values as unspecified and assigns zero
objective/gradient contribution to them.

## Lp rotation

The iterative reweighted least-squares structure is retained. The convergence
check compares successive rotated loading matrices directly; this is equivalent
to the intended upstream check and avoids repeatedly reconstructing the prior
loading matrix through a matrix inverse.

## Diagnostics

`calc_fitstats` includes a self-contained regularized incomplete-gamma routine
and a Poisson-mixture noncentral chi-square CDF for RMSEA confidence limits.
This removes the R `pchisq`/`uniroot` dependency.

## Validation

`test_criteria` checks analytical gradients against centered finite differences
for the smooth criteria, including varimax, quartimin, Crawford-Ferguson,
geomin, entropy, oblimin, infomax, McCammon, Tandem I/II, oblimax, Bentler,
binormamin, varimin, bifactor and bi-geomin.

`test_rotation` checks:

- orthogonal convergence and orthonormal transformation matrices;
- oblique reconstruction and unit-diagonal factor correlations;
- objective improvement from the identity start;
- Haar-style random-start orthogonality;
- Kaiser normalization;
- EIV and echelon identifying constraints;
- Lp execution;
- exact-zero SRMR for a model-implied correlation matrix;
- simplicity diagnostics.
