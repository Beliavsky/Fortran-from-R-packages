# evd-fortran

Modern Fortran/FPM translation of the computational core of the R package
**evd 2.3-7.1** by Alec Stephenson and contributors.

The upstream package is an implementation of univariate, bivariate and
multivariate extreme-value distributions, simulation, dependence functions,
threshold methods, stochastic-process utilities and likelihood fitting.
This port deliberately omits plotting, S3/R object presentation, formula
parsing, profile-plot infrastructure and other user-interface code.

## Implemented computational areas

- Univariate GEV, Gumbel, Frechet, reversed/negative Weibull, GPD and
  two-component Gumbel-max (`gumbelx`) density/CDF/quantile/RNG routines.
- Probability transformations used by `extreme` and `order` distributions.
- R-compatible `mtransform` exponential-measure transformations.
- Nine bivariate extreme-value families: logistic, asymmetric logistic,
  Huesler-Reiss, negative logistic, asymmetric negative logistic, bilogistic,
  negative bilogistic, Coles-Tawn and asymmetric mixed.
- Bivariate probability, density, Pickands dependence (`A`), angular density
  (`h`), conditional-copula and simulation routines.
- Multivariate logistic and asymmetric-logistic probability, density,
  dependence and simulation routines.
- Censored and Poisson-process bivariate POT likelihoods.
- Nonparametric bivariate CFG, Pickands, Tiago-de-Oliveira and POT dependence
  estimators, multivariate Pickands estimator, and exponential-scale contour
  construction.
- Clustering, extremal-index estimation, MARMA/MAR/MMA processes and `evmc`
  simulation.
- Stationary maximum-likelihood fitting for GEV, GPD, point-process,
  Gumbel-max and the nine bivariate maxima families.

The supplied `r_mod.f90` is reused for applicable R-compatible helpers and
optimization. No duplicate normal/beta/RNG/optimizer implementations were
added.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

BLAS and LAPACK are linked because the supplied `r_mod.f90` includes helpers
that use them.

The translation was also compiled directly with GNU Fortran 14.2 using
Fortran 2018 and `-Werror=implicit-interface`.

## Main modules

- `evd_univariate`
- `evd_transform`
- `evd_bivariate`
- `evd_multivariate`
- `evd_simulation`
- `evd_process`
- `evd_nonparametric`
- `evd_bvpot`
- `evd_fit`
- `evd` (convenience umbrella module)

## License and attribution

The evd-derived translation is distributed under the same **GPL-3** terms as
the upstream package. The separately supplied `src/r_mod.f90` remains under
its MIT license. See `LICENSES/`, `NOTICE`, and `upstream/`.

Please cite the original package when its algorithms are used:

A. G. Stephenson (2002), *evd: Extreme Value Distributions*, R News 2(2),
31-32.
