# Porting notes

Upstream `fastmatrix` 0.6-6 is GPL-3 and contains R, C, fixed-form Fortran and BLAS/LAPACK-facing code. This port rewrites the computational layer in modern free-form Fortran and does not require R headers.

## Algorithm choices

- LU/LDL, Gaussian elimination, Jacobi eigenanalysis, CG, Gauss-Seidel, power iteration, Sherman-Morrison, sweep, Householder, Floyd-Warshall and de Casteljau are direct self-contained modern implementations.
- Symmetric matrix square root, whitening, exponential, logarithm and power are evaluated through the Jacobi spectral decomposition.
- OLS solves the normal equations in the standalone implementation. `ols_fit_cg` provides the iterative normal-equation path corresponding to the upstream CG option.
- `rmnorm` uses a Cholesky factor and standard-normal Box-Muller draws.
- Chi probabilities are evaluated with a native regularized incomplete-gamma routine and quantiles by monotone bisection.

## Upstream provenance

The original `src/`, `DESCRIPTION`, and `NAMESPACE` are retained under `orig/` for license/provenance and parity work. They are not compiled by FPM.

## v0.2.0 compatibility additions

The lower-priority numerical parity targets from v0.1.0 are now implemented:

- real nonsymmetric Schur decomposition through LAPACK `DGEES`, including real 2x2 Schur blocks for complex conjugate eigenpairs;
- thin SVD through LAPACK `DGESVD`;
- QR least squares through LAPACK `DGELS`;
- SVD least squares with numerical-rank handling and covariance reconstruction;
- Mardia multivariate skewness/kurtosis coefficients and normality statistics/p-values;
- Harris variance-homogeneity tests in `Wald`, `log`, `robust`, and `log-robust` forms;
- callback-based Parlett recurrence matching upstream `matrix.fun` for upper-triangular matrices.

Unlike the v0.1.0 core, these compatibility routines intentionally link LAPACK/BLAS, as the upstream package itself does. The R print/S3 wrappers are represented by Fortran result types rather than emulated.

## Remaining nonessential parity targets

What remains is mostly R-facing infrastructure or aliases: formula/model-frame construction, print/summary/logLik/deviance methods, object-class predicates, and small `*.info`/extraction helpers whose numerical content is already available through the underlying Fortran arrays. The upstream `matrix.fun` itself is triangular-only; the earlier v0.1.0 note calling it nonsymmetric was inaccurate.

