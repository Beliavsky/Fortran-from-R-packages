# joker-fortran

Modern Fortran translation of the numerical core of the R package `joker`
(version 0.14.2).

## Scope

The port provides:

- d/p/q/r routines for Bernoulli, beta, binomial, Cauchy, chi-square,
  exponential, F, gamma, geometric, Laplace, log-normal, negative-binomial,
  normal, Poisson, Student-t, uniform, and Weibull distributions.
- categorical, multinomial, Dirichlet, and cumulative multivariate-gamma
  density/random/likelihood routines.
- the exported log-likelihood families from `joker`.
- closed-form MLE/ME estimators where the R package has them.
- beta, gamma, Weibull, Dirichlet and multivariate-gamma numerical MLEs.
- beta/gamma/Dirichlet/multigamma SAME or moment estimators used upstream.
- the complete upstream asymptotic covariance surface: MLE/ME formulas for
  all families that define them, plus SAME formulas for beta, gamma,
  Dirichlet, and multivariate-gamma.
- native incomplete beta/gamma, normal quantile, digamma, trigamma and
  inverse-digamma support.

The R S4 class system is intentionally represented as explicit Fortran
procedures and small derived result types rather than emulated.

## Build

```sh
fpm test
```

or compile directly with a Fortran 2018 compiler. Source lines conform to the
standard free-form 132-column limit; no compiler-specific line-length flag is
required.

## Sample-major convention

For multivariate observations, Fortran arrays use one observation per row:
`x(nobs, ndim)`.  This differs from a few upstream R helpers that operate on
category-by-sample matrices.

## License

GPL-3.0-or-later, matching upstream `License: GPL (>= 3)`.

## Asymptotic covariance API

Version 0.2.0 adds `joker_asymptotics`, including Fisher-information helpers
and estimator asymptotic covariances. The upstream coverage is complete for
the families that define `avar_*` methods:

- Bernoulli, binomial, categorical, chi-square, exponential, geometric,
  negative-binomial, Poisson: MLE and ME.
- Beta and gamma: MLE, ME, and SAME.
- Cauchy: MLE.
- Laplace, log-normal, and normal: MLE and ME.
- Dirichlet and multivariate-gamma: MLE, ME, and SAME.
- Multinomial: MLE and ME.

Fisher-F, Student-t, uniform, and Weibull have no `avar_mle`, `avar_me`, or
`avar_same` methods in upstream joker 0.14.2, so no non-upstream covariance
formula is invented for them.
