# Compatibility notes

The Fortran API maps all 12 exported computational R functions, but package
status is **substantial**, not complete R compatibility.

## Deliberate interface adaptations

- R `Matrix` classes are represented by the public `sparse_matrix` COO type.
- Every solver accepts either a vector right-hand side or a dense rank-2 `B`;
  matrix right-hand sides are solved column by column, matching the upstream R
  wrapper strategy. The aggregate `solver_matrix_result` retains each column's
  full error history and flags. Sparse R `Matrix` objects used as `B` are not a
  distinct Fortran type; callers pass their numerical columns as `real(dp)`.
- The Fortran matrix-RHS BiCG overload deliberately corrects an upstream R
  wrapper defect: `lsolve.bicg()` requests `errors1` and `errors2`, although the
  C++ kernel returns only `errors`.
- The upstream default `xinit=NA` draws `rnorm()` values, and even replaces an
  already-converged initial point by another random vector. The Fortran API
  defaults to zero and uses a deterministic sine vector only in that unusual
  already-converged case, making tests and runs reproducible.
- Stationary methods follow upstream behavior by perturbing exact zero diagonal
  entries; the translation uses deterministic `1e-5` instead of a random
  normal perturbation divided by `1e5`.
- Sparse matrix-vector products stay sparse. Sparse user preconditioners are
  converted once to dense storage for the reusable portable LU solve. Sparse
  SSOR triangular factors and sparse overdetermined normal equations currently
  materialize dense matrices.

## Numerical adaptations

- Chebyshev can accept `eigmin` and `eigmax` explicitly. When omitted, spectral
  endpoints are estimated by power and inverse-power iterations so the package
  remains self-contained; upstream uses R's `eigen(solve(M,A))`.
- QMR reproduces the upstream recurrence and flag meanings. For a non-identity
  preconditioner it forms a no-pivot L/U split, because the Armadillo/Matrix LU
  object conventions (including permutations) do not map directly to the two
  triangular factors consumed by the upstream recurrence.
- Breakdown comparisons use machine-safe near-zero tests rather than exact
  floating-point equality.

Underdetermined systems remain unsupported, matching the upstream public R
wrappers. Overdetermined systems are transformed to normal equations as in R.
