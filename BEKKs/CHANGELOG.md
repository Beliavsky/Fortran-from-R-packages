# Changelog

## 1.4.7-fortran.1

- Translated full, diagonal, and scalar BEKK models.
- Added symmetric and asymmetric Gaussian QML filtering and estimation.
- Added BHHH optimization, numerical score/Hessian inference, OPG covariance, and QML sandwich covariance.
- Added deterministic and randomized starting-value searches.
- Added simulation, fixed innovations, stationarity checks, and unconditional covariance.
- Added multi-step forecasts and parameter-uncertainty bands.
- Added volatility impulse responses with numerical delta-method bands.
- Added marginal and portfolio VaR, Kupiec and Christoffersen tests, and rolling backtests.
- Added portmanteau diagnostics and Monte Carlo parameter-recovery evaluation.
- Added elimination, duplication, commutation, selection, cut, lag, square-root, and pseudoinverse utilities.
- Corrected the portmanteau p-value calculation.
- Removed the constant covariance intercept from VIRF differences, where it cancels analytically.
- Added FPM metadata, applications, examples, tests, reproducible GNU Fortran builds, licensing, and provenance documentation.
