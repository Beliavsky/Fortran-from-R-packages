# Changelog

## 0.1.0 - 2026-08-20

- Initial modern Fortran/FPM translation of ZeroOneDists 1.0.0.
- Added BER, BER2, UHLG, UMB, and UPHN d/p/q/r routines.
- Added scalar and vector-parameter random generation APIs.
- Added GAMLSS-style family links, scores, working Hessians, deviance,
  initialization, means, and validity helpers.
- Added generic design-matrix regression fitting with BFGS.
- Stabilized the BER2 uniform boundary at `mu=0.5, nu=1`.
- Added independent reference, inversion, derivative, RNG, and fitting tests.
