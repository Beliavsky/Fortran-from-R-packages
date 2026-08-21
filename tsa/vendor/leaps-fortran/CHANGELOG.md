# Changelog

## 0.1.0

- Initial FPM translation of the `leaps` 3.2 computational core.
- Ported the R-level subset-selection setup and result calculations to a
  Fortran module API.
- Integrated the supplied Alan Miller free-format Fortran 95 QR and
  subset-search implementations.
- Added exhaustive, backward, forward, and sequential-replacement searches.
- Added forced-variable masks, rank-deficiency handling, model metrics,
  coefficient recovery, covariance matrices, examples, and tests.
- Omitted plotting and R-specific object/formula infrastructure.
