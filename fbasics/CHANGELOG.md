# Changelog

## 0.3.0

- Added VAR(p) moment-condition prewhitening with no-intercept least-squares estimation.
- Added exact long-run recoloring using `(I - A_1 - ... - A_p)^{-1}`.
- Added bounded conditional-sum-of-squares ARMA(1,1) fitting.
- Added the Andrews ARMA(1,1) plug-in bandwidth formulas for all five HAC kernels.
- Added automatic `newey_west`, `andrews_ar1`, and `andrews_arma11` bandwidth methods.
- Integrated prewhitening and automatic bandwidth selection into identity, two-step, iterated, and CUE GMM covariance calculations.
- Added bandwidth and prewhitening diagnostics to `gmm_result`.
- Added a fifth strict test suite covering ARMA fitting, formula reconstruction, VAR(1)/VAR(2) prewhitening, exact recoloring, GMM integration, PSD checks, and fallback behavior.

## 0.2.0

- Added Nolan-S1 stable density, CDF, quantile, RNG, ECF fitting, bounded MLE, Hessian, and covariance calculations.
- Added FMKL and five-parameter generalized-lambda distributions and robust, MLE, MPS, GOF, and histogram fitting paths.
- Added GH, HYP, GHT, SGH, SNIG, and SGHT fit wrappers and robust quantile-moment summaries.
- Added penalized B-spline smoothing-density estimation with density, CDF, quantile, and RNG procedures.
- Added identity, two-step, iterated, and CUE GMM.
- Added EL, ET, CUE, and ETEL generalized empirical likelihood.
- Added IID/HAC covariance calculations, five kernels including Quadratic Spectral, Newey-West and Andrews AR(1) bandwidths, J tests, and linear restriction tests.
- Added triangulated piecewise-linear interpolation and ordinary kriging with fitted exponential, Gaussian, or spherical covariance models.
- Added Shapiro-Wilk, Ansari-Bradley, Mood, Bartlett, Fligner-Killeen, Wilcoxon, Kruskal-Wallis, one-sample KS/Pearson, and adjusted Jarque-Bera calculations.
- Expanded the CSV application with FMKL, FM5, GH-family, stable, and spline-density modes.
- Added a fourth strict numerical regression suite covering all new algorithms.

## 0.1.0

- Added modern Fortran matrix, lag, PDL, positive-definite, and vectorization utilities.
- Added descriptive, row-wise, robust-quantile, L-moment, covariance, and correlation statistics.
- Added reproducible LCG and common continuous random generators.
- Added Normal, Student-t, NIG, GH, standardized GH, hyperbolic, GHT, and RS-GLD numerical routines.
- Added Normal, Student-t, NIG, and GLD fitting workflows.
- Added interpolation, hypothesis-test, and maximum-drawdown calculations.
- Added demo, CSV analyzer, strict builds, numerical tests, provenance, API mapping, and GPL header enforcement.
