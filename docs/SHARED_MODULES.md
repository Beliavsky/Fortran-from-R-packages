# Shared numerical modules

The translations contain many independent implementations of common R-like
operations. `rfortran-core` is a small MIT-licensed FPM library intended to
centralize only behavior that can be specified and tested consistently across
packages.

## Design rules

- Keep the core independent of BLAS, LAPACK, operating-system APIs and
  package-specific status types.
- Put each coherent area in its own Fortran module. Clients should use
  `only:` imports from the defining module.
- Preserve `r_mod` as a migration facade, not as the primary API for new code.
- Specify missing-value, empty-input, denominator and endpoint behavior.
- Return status values for invalid numerical inputs instead of stopping.
- Require deterministic unit tests and comparisons against R for every shared
  behavior used by a translation.
- Migrate one helper family at a time and retain package-level compatibility
  wrappers when an existing public API has different edge-case behavior.

## Initial module graph

```text
r_kinds       r_status
    |             |
    +------ r_missing
    |             |
    +------ r_descriptive
                  |
             r_time_series

r_vectors  (lagged and repeated vector/matrix differencing)

r_distributions  (normal and central t/chi-square/F densities, CDFs, and quantiles)
r_ordering       (stable order indices and value sorting)
r_sorting        (average ranks and compatibility exports)
r_quantiles      (type-7 median and unweighted/weighted quantile estimators)
r_robust         (median absolute deviation)
r_stability      (log-sum-exp and log-mean-exp)
r_transforms     (log1p, expm1, log1mexp, softplus, logistic, logit)
r_special        (digamma, trigamma, log-beta, regularized gamma/beta, integer combinatorics)
r_rolling        (valid-only and right-aligned trailing reductions and moments)

r_mod  (compatibility facade over the modules above)

rfortran-linalg  (optional checked BLAS/LAPACK-backed linear algebra)
```

## Current migrations

- `FinTS` uses the shared missing-aware count, mean and variance routines and
  the shared univariate and multivariate autocovariance/autocorrelation
  implementations.
- `fracdiff` uses the shared mean implementation in fractional differencing
  and semiparametric estimators.
- `vrtest` retains its existing helper API while delegating means, variances,
  and autocorrelations to the shared modules.
- `tseries` delegates means, variances, autocovariances, and price
  differencing while retaining its established edge-case behavior.
- `rugarch` delegates its normal density and CDF while retaining its public
  distribution API and local quantile implementation.
- `bayesgarch`, `bzinb`, and `gkwdist` delegate compatible digamma,
  trigamma, and log-beta implementations while retaining their public names
  and package-specific invalid-input behavior.
- `ACDm`, `mfGARCH`, and `frapo` delegate complete-window trailing means while
  retaining their valid-only, NaN-padded, and optional-trimming interfaces.
- `argus`, `DiscreteWeibull`, `bayesm`, and `gkwdist` delegate stable
  elementary transforms. Their compatibility wrappers preserve public names
  and, for `DiscreteWeibull`, its historical invalid-domain sentinel.
- `boot`, `fattailsr`, and `bayesm` delegate stable logistic transforms while
  retaining their endpoint, finite-sentinel, or probability-clamping APIs.
- `fportfolio` delegates rolling sample standard deviation while preserving
  its leading-zero result, and PortfolioTesteR uses shared rolling mean and
  standard deviation for complete-window Bollinger bands.
- `performanceanalytics` delegates valid-only rolling correlation while its
  wrapper preserves the package's zero result for constant finite windows.
  PortfolioTesteR delegates right-aligned rolling correlation while preserving
  return alignment, pairwise finite filtering, and undefined-variance NaNs.
- PortfolioTesteR's stochastic-price and stochastic-RSI indicators delegate
  rolling minima and maxima while retaining their complete-finite-window
  requirement. `PINstimation` delegates the VPIN imbalance-window total to the
  linear-time shared rolling sum.
- `corpcor` delegates weighted means and reliability-unbiased variances while
  preserving its normalized-weight validation and package status codes.
  `glmnet` delegates its normal positive-weight mean and population-variance
  path while retaining legacy fallbacks and caller-supplied centers.
- `mixtools` and `GB2` delegate inverse weighted-ECDF quantiles while retaining
  their public names and package-specific status or error behavior.
- `fitdistrplus`, `quarks`, `isotone`, and `survey` retain their public
  quantile APIs while delegating their distinct interpolation conventions to
  `r_quantiles`. The shared module names each convention explicitly rather than
  treating weighted quantiles as a single estimator. Survey's twelve selectable
  rules and their constants now also live in the shared module.
- `ACDm`, `boot`, `bayesm`, `extraDistr`, `performanceanalytics`, and `survey`
  delegate inverse-normal calculations to `r_qnorm`, with thin wrappers where
  finite endpoint sentinels or probability clamping are part of the existing
  package API.
- `ACDm`, `fitdistrplus`, `extraDistr`, and `survey` delegate compatible
  regularized gamma or beta calculations while preserving package validation
  and endpoint behavior. GB2 retains its local regularized-beta kernel because
  a numerically tiny replacement difference perturbed a sensitive optimizer
  regression, but it delegates its compatible log-beta calculation.
- `bayesm`, `corpcor`, `mixtools`, and `performanceanalytics` delegate medians
  to `r_quantiles`; PerformanceAnalytics also delegates its public sort wrapper
  to the stable ordering module. ACDm delegates its empirical type-7 quantile.
- `bayesm`, `extraDistr`, `chyper`, and `BiasedUrn` delegate integer
  log-factorial or log-binomial-coefficient calculations while retaining
  their established invalid-input sentinels.
- `bayesm`, `bzinb`, `extraDistr`, and `mixtools` delegate log-sum-exp
  reductions to `r_stability`. Compatibility wrappers preserve bayesm's
  finite empty-vector sentinel and bzinb's historical absent-component
  sentinel.
- `waveslim` delegates finite-filtered means, variances, medians, MADs, and
  type-7 quantiles, as well as normal CDF/quantile and regularized gamma.
  Its lag-dependent autocovariance and cross-correlation remain local because
  their denominator policies differ from the R-style core time-series API.
- The central Student-t, chi-square, and F distribution families now share
  density, lower/upper-tail probability, log-probability, and quantile APIs.
  Waveslim's chi-square compatibility wrapper delegates to this public layer.
- `spantest` delegates normal probabilities and quantiles, Student-t
  probabilities, and F upper tails while retaining its invalid-input and
  finite-endpoint sentinels. `rrcov` delegates normal, chi-square, F, and
  regularized gamma/beta calculations while retaining its probability clamp
  and public helper names. `survey` delegates central t quantiles and F and
  chi-square survival probabilities while retaining its validation behavior.
- `vares` delegates its central Student-t and F density, probability, and
  quantile kernels. Its public distribution procedures retain their original
  names, optional defaults, log/tail flags, elemental interfaces, and finite
  Student-t endpoint sentinels; derived asymmetric-t and half-t families use
  the shared kernels through the existing compatibility layer.
- The complete test suites for the first five packages and for `bayesgarch`,
  `bzinb`, `gkwdist`, `ACDm`, `mfGARCH`, `argus`, `DiscreteWeibull`, and
  `bayesm`, `boot`, `fattailsr`, PortfolioTesteR, `performanceanalytics`,
  `PINstimation`, `corpcor`, `glmnet`, `mixtools`, `GB2`, `extraDistr`,
  `chyper`, `BiasedUrn`, `fitdistrplus`, `quarks`, `isotone`, and `survey` pass
  after migration. The `spantest`, `rrcov`, and `vares` suites also pass after
  their probability helpers were migrated.
  FRAPO's
  library and test executable compile, but its test process currently exits
  with a Windows stack-overflow code on this machine after linking against the
  available Octave BLAS/LAPACK libraries.
- `fportfolio` also compiles fully, but its test executables encounter the
  same Windows stack-overflow exit with those BLAS/LAPACK libraries.
- The separate `rfortran-linalg` package provides checked square solves,
  symmetric eigendecomposition, Cholesky factors, and SPD inverses with log
  determinants through a pinned pure-Fortran LAPACK backend. `CEoptim`,
  `cmaes`, and `cccp` delegate their compatible operations to this layer;
  `cccp` uses the intrinsic `norm2` for its Euclidean-norm compatibility API.
  All 15 package tests pass after these migrations, as do the direct shared
  linear-algebra tests.
- Their deterministic comparisons against R pass 7/7 and 9/9 cases,
  respectively, for `FinTS` and `fracdiff`; `rugarch` passes 25/25 cases.
- The direct `rfortran-core` comparison against R references passes all 54/54
  special-function, transform, descriptive, quantile, and rolling-statistic
  cases, including stable ordering, inverse-normal tails, median/MAD, and
  regularized gamma/beta, central t/chi-square/F distributions, and log-sum-exp
  reductions. The integrated comparison runner covers 151 cases.

## Planned layers

Possible later distribution additions include noncentral families, but those
require separate algorithms and substantially broader tail validation.
Optimization, higher-level models, and data-frame or file APIs should remain
separate optional libraries so lightweight translations do not inherit
unnecessary dependencies. BLAS/LAPACK operations now begin this separation in
`rfortran-linalg`; its API should expand only when multiple translations share
the same shape, failure, and numerical semantics.

Before deleting an old helper, check public imports, internal calls,
transpiler-generated calls, package tests and comparison cases. Procedures
unused by the current translations may still be useful to the transpiler and
should normally be made private or moved to an optional layer before removal.

## Dependency layout

During development, translated packages use a relative FPM path dependency on
`../rfortran-core`. Once the API stabilizes, the core should be published as a
separately tagged repository and package manifests should pin a Git revision.
That will allow an individual translation to be downloaded and built without
the full aggregate repository.

The optional `rfortran-linalg` package follows the same local-path arrangement
and pins its `fortran-lapack` dependency to an exact Git revision. A cold FPM
build compiles that backend separately for each independent translated package;
later builds reuse the package's local build products.
