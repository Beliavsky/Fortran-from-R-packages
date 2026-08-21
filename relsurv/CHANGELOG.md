# Changelog

## 0.2.0

- Added `rs.br` Brownian-bridge and Cramer-von Mises diagnostics.
- Added `rs.zph` time transforms and covariance scalings.
- Added `residuals.rsadd`-style Schoenfeld/partial-residual covariance kernels,
  including tied-event averaging.
- Added exact `epa` adaptive bandwidths, boundary kernels, and left-predictable
  smoothing.
- Added the `rsadd` EM branch with unknown-cause weights, weighted Cox M-steps,
  smoothed baseline hazards, automatic Ederer-II bandwidth selection, and
  missing-information covariance correction.
- Added grouped binomial and Poisson `rsadd` GLM branches.
- Added default, YL2013, and YL2017 years-lost calculations with Greenwood area
  variance and bootstrap-variance aggregation.
- Added `transrate.hld` and `transrate.hmd` parsers.
- Expanded the strict suite to 13 tests plus the example.

## 0.1.0

- Initial modern free-format Fortran/FPM translation of the computational core
  of `relsurv` 2.3-3.
- Translated population rate-table stepping and expected-survival machinery.
- Added Pohar-Perme, Ederer I/II, and Hakulinen net-survival estimators.
- Added crude-mortality, grouped net-survival comparison, and `nessie` routines.
- Translated Aalen additive-hazard numerical kernels.
- Added additive, multiplicative, and transformed-time regression paths.
- Added counting-process splitting, rate-table utilities, smoothing, and
  years-lost integration primitives.
- Integrated the supplied Fortran `survival` Cox implementation.
- Added strict regression tests and an FPM example.
