# MNB-fortran

Modern free-format Fortran/FPM translation of the computational layer of R package
`MNB` 1.2.0, *Diagnostic Tools for a Multivariate Negative Binomial Regression Model*.

The upstream package models clustered counts through a Gamma/GLG random intercept,
which produces the multivariate negative-binomial likelihood after integrating out
the random effect.

## Build

```text
fpm build
fpm test
fpm run --example demo_mnb
```

The package can also be compiled directly with GNU Fortran; see `scripts/validate.sh`
and `scripts/validate.bat`.

## Main API

```fortran
use mnb

type(mnb_fit_result) :: fit
fit = fit_mnb(start, y, x, n, mi, offset)
```

Inputs are explicit numerical arrays rather than R formulas/model frames:

- `y(n*mi)` contains clustered counts in subject-major order.
- `x(n*mi,p)` is the already constructed design matrix, including an intercept if desired.
- `n` is the number of subjects/clusters.
- `mi` is the common number of observations per subject, matching the upstream package assumption.
- `offset` is optional and is supplied on the linear-predictor scale.
- parameter vectors are ordered as `phi, beta(1:p)`.

Implemented computational functionality includes:

- integrated multivariate negative-binomial log likelihood;
- BFGS maximum-likelihood fitting, numerical Hessian, covariance, standard errors,
  z statistics and two-sided normal p-values;
- simulation from the Gamma-Poisson representation used by upstream `rMNB`;
- randomized quantile residuals;
- weighted, standardized weighted, Pearson, standardized Pearson and cluster
  deviance residuals;
- cluster-deletion generalized Cook distances and likelihood displacement;
- local influence under subject case weights, observation case weights,
  covariate perturbation and dispersion perturbation;
- Atkinson-style Monte Carlo residual envelopes.

Plotting and R formula/S3/data-frame infrastructure are deliberately omitted.

## Source-compatible quirks

The Fortran implementation preserves two numerical behaviors of MNB 1.2.0 that
appear unusual but affect exact source parity:

1. `re.MNB` constructs `W` with a global outer product of the fitted means rather
   than a block-diagonal cluster covariance. `residuals_mnb` reproduces that formula.
2. `cova.pertu` constructs a perturbed `X.new` but evaluates the likelihood with the
   original `X`. Consequently the source-compatible covariate local-influence block
   is numerically zero. The Fortran port preserves this behavior rather than silently
   changing the diagnostic.

The local-influence routine also preserves upstream selection of the *last* eigenvector
returned by R `eigen()`, even though the documentation calls it the maximum-curvature
direction.

## Dependencies

The compiled library is self-contained. The only numerical uses of upstream package
imports are narrow:

- `flexsurv::rgengamma(mu=0, sigma=1/sqrt(phi), Q=1/sqrt(phi))` simplifies exactly
  to a Gamma random effect with shape `phi` and scale `1/phi`, implemented natively.
- `numDeriv` Hessians are represented by native central finite differences.

The available Fortran translations of `flexsurv` and `numDeriv` are retained under
`vendor/` for provenance and future parity work.
