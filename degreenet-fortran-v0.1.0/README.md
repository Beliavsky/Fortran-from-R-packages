# degreenet-fortran

Modern Fortran 2018/FPM translation of the computational core of the R package
`degreenet` 1.3-7 (Mark S. Handcock / statnet).

**Required upstream attribution:** Based on 'statnet' project software
(statnet.org). For license and citation information see statnet.org/attribution.

## Scope

The translation covers the numerical/statistical code used by the package:

- Riemann zeta via the upstream Euler-Maclaurin construction.
- Yule, discrete Pareto/Zipf, tapered power, Waring, discrete q-exponential,
  generalized hypergeometric degree (GHDI), negative-binomial and
  exponential-power degree laws.
- Conway-Maxwell-Poisson normalization, natural-parameter PMF, random
  generation support, mean/SD conversion, and mean/SD parameterization.
- Poisson-lognormal PMF and simulation.
- Geometric-Yule, geometric-discrete-Pareto, negative-binomial-Yule,
  negative-binomial-Waring and geometric-Waring stopped-process laws.
- Untruncated/right-truncated likelihood evaluation for all model families.
- A common bounded Nelder-Mead MLE engine with numerical Hessian,
  covariance matrix and standard errors, replacing the repeated R `optim()`
  wrappers.
- Grouped-degree probabilities corresponding to the package's codes 5--10.
- Rounded-degree observation bins used by the `llr*` routines.
- Modified Anderson-Darling statistic, penalized Hellinger distance and the
  package concentration-index calculation.
- Parametric model simulation and bootstrap fitting.
- Simple graph realization from a degree sequence and Yule-degree graph
  generation without `igraph`/`network`.

Plotting, R S3/list formatting, file-saving side effects, package load hooks,
and R-specific `igraph`/`network` objects are intentionally omitted.

## Build

```text
fpm build
fpm test
fpm run --example fit_degree_models
```

The code has no external Fortran dependencies.

## Parameterization notes

The low-level `dwar()` translation follows the exported R density and therefore
uses Waring natural parameters `(rho, a)`.  R likelihood wrappers such as
`awarmle()` transform their second user parameter (`prob. new`) to `a`; use
`waring_prob_to_natural()` before fitting when reproducing that interface.

`MODEL_CMP` uses CMP natural parameters `(lambda, nu)`.  Use
`cmp_mutonatural()` / `cmp_moments()` for the R package's mean/SD interface.

The common model constants and `model_pmf()`, `loglik_model()` and
`fit_degree_model()` APIs replace the large family of nearly identical R
`ll*` / `a*mle` wrappers.

## License

GPL-3.0-or-later with the additional statnet attribution requirements in the
upstream `LICENSE`, retained verbatim in this distribution.
