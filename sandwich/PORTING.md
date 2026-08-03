# Porting notes

## Computational representation

The R package is organized around `estfun()` and `bread()` S3 generics. The
Fortran translation moves this interface boundary outward:

- an estimating-function matrix is `scores(n, k)`;
- a bread matrix is `bread(k, k)`;
- cluster and time identifiers are explicit integer arrays;
- covariance estimators return allocatable real matrices.

This representation preserves the numerical algorithms while avoiding a
Fortran reimplementation of R formulas, environments, model frames, and class
dispatch. The supplied `ols_model` is the built-in adapter for weighted linear
least squares. Other model implementations can produce scores and bread using
their own derivatives.

## Scaling convention

The bread is scaled as `n * inverse(information)`. Meat procedures return
observation-averaged cross-products. The final covariance is
`bread * meat * bread / n`, matching the upstream package's convention.

## Algorithm mapping

- `sandwich.R`: `meat`, `sandwich_covariance`, and `vcov_opg`.
- `bread.R` and `estfun.R`: generic matrix interface plus OLS adapter.
- `vcovHC.R`: `hc_weights`, `meat_hc`, and `vcov_hc`.
- `vcovHAC.R`: kernel functions, HAC meat, VAR prewhitening, Andrews and
  Newey-West bandwidths, and Lumley weights.
- `lrvar.R`: `long_run_variance` through an intercept-only score calculation.
- `auxiliary.R`: PAVA, autocorrelation, and isotonic ACF.
- `vcovCL.R`: one-way and multiway clustered covariance.
- `vcovPL.R`: panel longitudinal/Driscoll-Kraay covariance.
- `vcovPC.R`: panel-corrected covariance.
- `vcovBS.R`, `vcovBS.lm.R`, and `vcovJK.R`: general replicate covariance and
  OLS cluster/bootstrap implementations.

## Numerical choices

The project is self-contained. Linear systems use scaled partial-pivoting
Gaussian elimination. Symmetric eigendecompositions use a Jacobi method, which
supports the cluster HC2 matrix powers and PSD repair without LAPACK.
Andrews ARMA(1,1) bandwidth fitting uses a bounded conditional-sum-of-squares
Nelder-Mead search rather than R's `arima()` maximum-likelihood implementation.
The AR(1) path and kernel/bandwidth formulas follow the upstream code directly.

## Deliberate omissions

The following are R infrastructure rather than portable computational kernels:

- S3 method dispatch and adapters for `glm`, `coxph`, `nls`, `polr`, `mlogit`,
  `survreg`, `rlm`, hurdle, zero-inflated, and other external model classes;
- formula evaluation, model frames, `na.action`, `zoo` indexing, names, and
  attributes;
- ordering by an R formula and data environment;
- R callback functions for user-defined omega, kernel, bootstrap, and refit
  methods;
- parallel `apply` backends;
- multiresponse `mlm` object formatting.

These omissions do not prevent use with those statistical models: a Fortran
model can call the generic robust covariance routines after supplying its score
and bread matrices.
