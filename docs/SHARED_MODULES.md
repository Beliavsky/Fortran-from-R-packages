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
- All four package test suites pass after migration.
- Their deterministic comparisons against R pass 7/7 and 9/9 cases,
  respectively.

## Planned layers

Likely additions are rolling statistics, `r_special`, `r_distributions`, and
`r_sorting`. BLAS/LAPACK linear algebra, optimization, random-number
generation, higher-level models, and data-frame or file APIs should remain
separate optional libraries so lightweight translations do not inherit
unnecessary dependencies.

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
