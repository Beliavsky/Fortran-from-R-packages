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

r_distributions  (normal density and CDF)
r_sorting        (merge sorting, type-7 quantiles, average ranks)
r_stability      (log-sum-exp and log-mean-exp)
r_special        (positive-domain digamma, trigamma, log-beta)
r_rolling        (valid-only and right-aligned trailing means)

r_mod  (compatibility facade over the modules above)
```

The first time-series pilot deliberately starts small:

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
- The complete test suites for the first five packages and for `bayesgarch`,
  `bzinb`, `gkwdist`, `ACDm`, and `mfGARCH` pass after migration. FRAPO's
  library and test executable compile, but its test process currently exits
  with a Windows stack-overflow code on this machine after linking against the
  available Octave BLAS/LAPACK libraries.
- Their deterministic comparisons against R pass 7/7 and 9/9 cases,
  respectively, for `FinTS` and `fracdiff`; `rugarch` passes 25/25 cases.
- The direct `rfortran-core` comparison against base R passes 5/5 digamma,
  trigamma, log-beta, and rolling-mean cases. The integrated comparison runner
  passes all 102/102 current cases.

## Planned layers

Likely additions are rolling variance and standard deviation, stable elementary
transforms such as `log1p` and `log1mexp`, and special functions only after
their tail and domain policies have been compared carefully. BLAS/LAPACK linear
algebra, optimization, random-number generation, higher-level models, and
data-frame or file APIs should remain separate optional libraries so
lightweight translations do not inherit unnecessary dependencies.

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
