# pan

Modern free-form Fortran translation of the computational code in the R package
`pan` 2.0 (CRAN publication 2026-06-30).

The upstream package implements multiple imputation for multivariate panel or
clustered data under a multivariate linear mixed model and maximum-likelihood
estimation for a univariate linear mixed model. This translation keeps the
numerical/statistical core and omits R-specific S3/list handling, `.Fortran`
marshalling, plotting, data-frame operations, and package-registration code.

## Implemented API

The public facade is the module `pan`:

```fortran
use pan, only : dp, pan_prior, pan_result, pan_mcmc
use pan, only : pan_bd_prior, pan_bd_result, pan_bd_mcmc
use pan, only : ecme_result, ecme_fit
```

### `pan_mcmc`

Implements the full multivariate mixed-model Gibbs sampler corresponding to
R `pan()`.

* arbitrary numbers of response variables, fixed effects, and random effects;
* full covariance across all vectorized random effects;
* inverse-Wishart updates for residual and random-effect covariance matrices;
* matrix-normal fixed-effect draws;
* Gaussian subject-specific random-effect draws;
* arbitrary component-wise missingness in the response matrix;
* conditional multivariate-normal imputation, including completely missing rows;
* deterministic package-local RNG;
* restart states equivalent in purpose to R `result$last`.

Missing response values are represented by IEEE NaNs. Predictors must be
complete. Subject labels must be sorted into contiguous blocks, as required by
the upstream R package.

### `pan_bd_mcmc`

Implements the `pan.bd()` model, where the random-effect covariance is
block-diagonal across response variables. Each response has its own `q x q`
inverse-Wishart-updated random-effect covariance block.

### `ecme_fit`

Fits the univariate Gaussian mixed model used by R `ecme()`:

```
y_i = X_i beta + Z_i b_i + e_i
b_i ~ N(0, Psi)
e_i ~ N(0, sigma2 V_i)
```

`V_i` is extracted from an optional maximum occasion covariance matrix
`vmax`; identity covariance is used when `vmax` is absent. With no `zcol`,
the fit is exact generalized least squares. With random effects, the
translation uses a Gaussian EM/ECME-target iteration with the same
maximum-likelihood objective and returns fixed effects, `sigma2`, `Psi`,
marginal covariance of beta, log likelihood, empirical Bayes random-effect
means, and conditional random-effect covariances.

## Numerical implementation

The upstream `pan.f` contains its own Cholesky, triangular solve, Wishart, and
random-number routines. The translation therefore remains self-contained and
does not vendor BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or another
translated R package. Dense SPD linear algebra is implemented directly in
modern Fortran and is adequate for the relatively small covariance dimensions
used by these algorithms.

All maintained real-valued Fortran code uses the single named kind `dp =
real64`. Explicit IEEE NaN checks use `ieee_is_nan`.

## Build

From the top-level `pan` directory:

```text
fpm build
fpm test
fpm run --example pan_example
```

No separately installed BLAS or LAPACK library is required.

## Example

`example/pan_example.f90` runs a small missing-data random-intercept Gibbs
sampler and an `ecme_fit` random-intercept model.

## Tests

The deterministic test suite covers:

* Park-Miller RNG reproducibility and inverse-Wishart SPD output;
* univariate and multivariate `pan_mcmc`;
* component-wise and completely missing rows;
* full random-effect covariance;
* `pan_bd_mcmc`;
* restart-state reproducibility;
* exact no-random-effects GLS;
* non-identity `vmax`;
* iterative random-intercept ML fitting;
* invalid subject ordering.

See `VERIFICATION.md` for the commands run in the translation environment.

## License and provenance

The upstream package is GPL-3. This translation is distributed under
GPL-3.0-only. See `LICENSE`, `NOTICE.md`, `PROVENANCE.md`, the retained
`UPSTREAM_DESCRIPTION`, and `UPSTREAM_CITATION`.
