# Translation coverage

## Directly represented computational behavior

The first Fortran release covers the principal numerical functionality exposed
by `gsl_nls()` and `gsl_nls_large()`:

| Upstream functionality | Fortran status |
|---|---|
| Levenberg-Marquardt trust region | Implemented |
| LM geodesic acceleration | Implemented |
| Dogleg | Implemented |
| Double dogleg | Implemented |
| 2-D subspace trust-region step | Implemented |
| Steihaug-Toint CG | Implemented |
| Matrix-free Jacobian operator path | Implemented |
| Parameter bounds | Implemented |
| Observation weights | Implemented |
| General positive-definite weight matrix | Implemented |
| Forward/central finite differences | Implemented |
| Directional second derivative | Implemented |
| Robust IRLS | Implemented |
| Huber/Barron/bisquare/Welsh/optimal/Hampel/GGW/LQQ | Implemented |
| Multi-start fitting | Implemented, with quasi-random difference noted below |
| Covariance/hat/Cook/log-likelihood helpers | Implemented numerically |
| Iteration traces | Implemented |

## Deliberate implementation differences

### No GSL runtime dependency

Upstream delegates low-level trust-region and linear algebra operations to GSL.
The Fortran package implements these operations directly. QR, Cholesky, a
Jacobi-eigen pseudoinverse, trust-region steps, and finite differences are all
self-contained. The algorithmic branches are modeled on the package/GSL
behavior, but floating-point trajectories and iteration counts need not be
bit-identical to GSL.

### Multi-start quasi-random generator

Upstream uses a GSL Sobol generator for parameter dimension below 41 and Halton
above that threshold. This release uses a Halton sequence for every dimension.
The Hickernell-Yuan-style retained-sample/local-search orchestration is retained,
but exact starting sequences can therefore differ from the R package for
`p < 41`.

### Sparse Matrix input

The R large-system front end accepts `Matrix` sparse Jacobians. The Fortran
translation instead provides `fit_nls_large_operator`, which is more general:
the caller supplies `J*v` and `J^T*v` without materializing either a dense or a
sparse matrix. A dedicated CSC object front end is not included.

### Parallel execution

Independent evaluations used by the R implementation may run through external
R/GSL execution machinery. This release is serial. The numerical objective and
trust-region logic are unchanged by that execution choice.

### R model and inference infrastructure

Not translated literally:

- formula parsing and environments;
- `deriv()`/`numericDeriv()` expression handling;
- `Matrix` S4 object conversion;
- S3 `nls`-compatible model objects;
- printing, summaries, ANOVA, prediction data-frame handling;
- profile likelihood and profile-based confidence intervals;
- the exported MINPACK/GSL test-problem catalog (`nls_test_list` and
  `nls_test_problem`) as R objects.

The original test catalog and source remain in `original/gslnls-master/`.

## Robust-loss fidelity

The tuning constants and psi/dpsi formulas are translated from
`src/nls_irls.c`. IRLS uses the upstream robust scale constant
`1.482602218505602 * median(abs(residual))`, normalizes robust weights to mean
one, and combines them with user-supplied observation weights.

## General weight matrices

The R source passes the lower factor produced by `t(chol(W))`; the translated
routine follows that lower-factor convention when transforming residuals and
Jacobians.

## Callback portability

All direct user callback invocations in library code occur in module-level
procedures with explicit `procedure(...)` dummy interfaces. In particular, the
matrix-free Steihaug routine is module-level rather than an internal procedure.
This avoids a gfortran portability failure seen with host-associated callback
procedures when `-Werror=implicit-interface` is enabled.
