# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational-core port of Rfast2 0.1.5.6.
- Integrated supplied `Rfast-fortran-v0.3.0` as a path dependency.
- Added Rfast2 array/group/quantile/trim/metric kernels.
- Added portable PCG32-based random generation and sampling.
- Added univariate/statistical tests, Kaplan-Meier, circular correlations,
  meta-analysis, silhouette, permutation/bootstrap and energy-test kernels.
- Added Rfast2-specific MLEs and column-wise MLE routines.
- Added logistic, Poisson, gamma, Weibull, multinomial, binomial, ZTP,
  tobit and related regression/scan interfaces.
- Added constrained/robust/heteroscedastic least squares.
- Added PCA/PCR, Mahalanobis depth, leverage, item discrimination/difficulty,
  and covariance-distance calculations.
- Corrected documented upstream RNG, sampling, Poisson-rate, gamma-Poisson,
  and numerical-stability defects.
