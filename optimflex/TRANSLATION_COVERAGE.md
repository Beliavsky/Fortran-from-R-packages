# Translation coverage

Upstream: `optimflex` 0.1.8.

All functions exported by the R package are represented in the Fortran API:
`bfgs`, `dogleg`, `double_dogleg`, `fast_grad`, `fast_hess`, `fast_jac`,
`gauss_newton`, `is_pd_fast`, `l_bfgs_b`, `levenberg_marquardt`,
`modified_newton`, and `newton_raphson`.

## Deliberate representation differences

- R closures and `...` are replaced by typed Fortran procedure callbacks.
- R control lists are replaced by `optim_control` plus method-specific default
  constructors.
- `numDeriv` is not required. The Richardson option is implemented directly in
  Fortran using Richardson-extrapolated centered differences. This is the same
  numerical technique, but it is not intended to reproduce every internal
  tuning constant of `numDeriv` bit-for-bit.
- `modified_newton` exposes the common `diff_method` field rather than R's
  separate `grad_diff`/`hess_diff` strings. Complex-step differentiation is not
  exposed because the public objective callback is real-valued.
- Dogleg/double-dogleg keep a dense symmetric BFGS curvature matrix. The R code
  stores the mathematically equivalent BFGS curvature through a Cholesky factor
  and rank-one update/downdate operations for additional numerical robustness.
- R exception handling and `NA` values map to explicit result status strings
  and IEEE finite-value checks.
- BFGS accepts `lower`/`upper` in the R signature but does not use them in the
  upstream implementation; the Fortran BFGS API therefore does not pretend to
  enforce box constraints. Use `l_bfgs_b`, dogleg, double-dogleg, or LM for
  bounds.

No plotting code exists in the package. R documentation/testthat infrastructure
is not translated.
