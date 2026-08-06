# Porting notes

## Scope

The port covers all exported computational routines and the native C++ kernels
in `spantest` 1.4-0. It omits Rcpp registration, roxygen/package hooks, R lists,
S3-style presentation, and testthat infrastructure. The upstream package has no
plotting routines.

## Linear algebra

R uses QR decompositions, `solve`, and RcppArmadillo. The Fortran port uses
partial-pivoting Gaussian elimination and normal-equation residualization. The
same estimands and restrictions are evaluated. Very ill-conditioned problems
may differ numerically; singular systems return an explicit status instead of
raising an R condition.

## Probability functions

Normal probabilities use the Fortran `erfc` intrinsic. Normal quantiles use the
Acklam rational approximation. Student-t and F probabilities use a
continued-fraction implementation of the regularized incomplete beta function.

## Randomized tests

`span_gl_a`, `span_gl_ad`, and randomized `span_as` runs use a local xorshift64
stream with Box-Muller normals. Seeds make runs deterministic within the
Fortran package without changing any process-global RNG state. Results are not
expected to match R's Mersenne-Twister draws bit for bit, but the sign-flip and
product-normal algorithms are the same.

The upstream C++ `gl_sim_stats` optimization is translated as a streaming
Fortran loop; it never materializes a `T x (N * nsim)` array.

## Subseries test

The upstream `cut(seq_len(T), k)` construction is represented by balanced,
contiguous integer folds. Batched Frisch-Waugh residualization, Student
subseries p-values, per-asset alpha/delta Cauchy merging, optional B-draw
merging, and final cross-sectional Cauchy merging are retained.

## Simulation

All twelve DGP presets are implemented. Toeplitz cross-sectional dependence is
introduced with an upper Cholesky factor. Student-t draws use normal and gamma
variates; the standardized Fernandez-Steel skew-t transformation follows the R
formula. The sequence of random draws differs from R, but seeded Fortran runs
are reproducible.
