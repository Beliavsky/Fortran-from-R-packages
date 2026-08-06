# Changelog

## 0.1.0 - 2026-08-04

- Initial modern Fortran 2018 computational translation of `tsgarch` 1.0.4.
- Added eight univariate volatility-model families.
- Added ten standardized innovation distributions through the vendored
  `tsdistributions` translation.
- Added filtering, bounded estimation, numerical inference, simulation,
  forecasting, PIT diagnostics, profiles, and rolling VaR backtests.
- Added variance targeting, variance regressors, initialization methods, exact
  IGARCH equality parameterization, and benchmark constants.
- Added FPM/Make builds, six tests, an example, API mapping, validation notes,
  and complete provenance archives.
- Retained the attached LGPL-3.0-or-later `nloptr` translation as unlinked
  provenance because it cannot be combined with GPL-2.0-only `tsgarch`.
