# Changelog

## 0.3.0

- Added `lmrob`-style S, SM/MM, and SMDM estimation chains.
- Added nonsingular, simple, exhaustive, and automatic S-subsampling modes.
- Added weighted, AVAR1-style, and sandwich regression covariance estimators.
- Added a tested `lmrob.lar` numerical implementation.
- Added `glmrob` Mqle binomial/Poisson and MT binomial/grouped-binomial estimators.
- Added `nlrob` MM, tau, constrained-M, and maximum-trimmed-likelihood estimators.
- Added continuous-data full adjusted outlyingness with sampled hyperplane directions.
- Added partitioned FAST-style MCD and LTS with full-data concentration refinement.
- Added analytical MCD consistency and finite-sample corrections.
- Added deterministic-MCD singular exact-fit hyperplane outputs.
- Added robust prediction intervals, Wald/deviance tests, robust R-squared, outlier statistics, and tolerance-ellipse coordinates.
- Added `full_rank_matrix` column reduction.
- Added Welsh, optimal, GGW, and LQQ score/loss/weight functions and derivatives.
- Added CLI modes `partlts`, `lar`, `lmrob`, `smdm`, `mqle-binomial`, `mqle-poisson`, and `mt`.
- Added bounds enforcement during nonlinear MM refinement.
- Added user-specified degrees of freedom to robust deviance comparisons.
- Expanded strict debug and optimized regression tests.

## 0.2.0

- Added `rankMM`/classical PCA, deterministic six-start MCD, advanced LTS search modes, and Bianco-Yohai logistic regression.

## 0.1.0

- Initial modern Fortran numerical translation of robust scales, score functions, medcouple, robust covariance analogues, LTS/MM regression, robust GLM, and nonlinear robust least squares.
