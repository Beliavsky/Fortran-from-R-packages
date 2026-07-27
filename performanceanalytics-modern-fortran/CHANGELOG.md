# Changelog

## 0.3.0

- Added exact analytical M2, M3, and M4 shrinkage workflows.
- Added exact sample MSE and target/sample covariance kernels for all original target families.
- Added the sixth-order k-statistic unbiased M3 MSE correction.
- Added independent/equal-marginal, observed one-factor, constant-correlation, Simaan, and central-symmetry targets with multiple-factor support.
- Added returned target vectors, A matrices, b vectors, constrained weights, sample moments, and corrected moments.
- Added bias-corrected structured independent/equal-marginal coskewness.
- Added `exact_moment_shrinkage` and a dedicated finite-sample regression suite.
- Corrected two apparent original C-kernel defects in distinct-index VM3 and all-distinct constant-correlation CM4 calculations.

## 0.2.0

- Added structured covariance, co-skewness, and co-kurtosis estimators.
- Added single- and multi-target M2/M3/M4 shrinkage.
- Added EWMA co-skewness and co-kurtosis.
- Added M3/M4 Moment Component Analysis using higher-order orthogonal iteration.
- Added constrained factor-model Nearest Comoment Estimation with PCA initialization and order weighting.
- Added GPD fitting, VaR, ES, standard errors, and confidence intervals.
- Added lognormal, Monte Carlo, and kernel portfolio risk engines.
- Added moving-block bootstrap risk standard errors.
- Added rolling, expanding, and conditional CAPM routines.
- Added capture ratios, trailing outperformance probabilities, and reusable summary statistics.
- Added `advanced_estimators` and a comprehensive advanced regression suite.

## 0.1.0

- Initial modern Fortran translation of returns, risk, drawdowns, performance ratios, CAPM, sample co-moments, portfolio moment risk, rolling statistics, cleaning, and applications.
