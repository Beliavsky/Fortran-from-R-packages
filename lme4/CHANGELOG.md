# Changelog

## 0.3.0

- Added `fit_lmm_pls`, a penalized-least-squares/Woodbury LMM solver that avoids the dense observation covariance matrix.
- Added diagonal, compound-symmetry, and AR(1) random-effect covariance structures alongside the existing unstructured form.
- Added `family_spec_t` callbacks and `fit_glmm_custom` for user-defined links, means, variances, likelihoods, and response validation.
- Added Gaussian-identity, binomial-probit, binomial-cloglog, and quasi-Poisson-log family specifications.
- Added multidimensional adaptive Gauss-Hermite quadrature for one grouped vector of correlated random coefficients.
- Added Gaussian nonlinear mixed-model fitting, prediction, and simulation with a user-provided nonlinear mean function.
- Added Wald and profile-likelihood fixed-effect confidence intervals.
- Added parametric bootstrap, percentile intervals, group-deletion influence diagnostics, Cook distances, and likelihood-ratio tests.
- Added `lmList`-style grouped weighted regression and prediction.
- Added an advanced example and expanded runtime tests for every new subsystem.

## 0.2.0

- Added Gamma-log and inverse-Gaussian-log GLMMs with fixed dispersion.
- Added fixed-size negative-binomial-log GLMMs.
- Added `fit_glmer_nb` bounded profiling of the negative-binomial size parameter.
- Added adaptive Gauss-Hermite fitting for one grouped scalar random-effect coefficient.
- Added Gamma, inverse-Gaussian, and negative-binomial simulation.
- Extended prediction and Pearson residuals to the new families.
- Added runtime tests for every new family, dispersion profiling, and binomial/Poisson adaptive quadrature.

## 0.1.0

- Added dense Gaussian LMM fitting by ML and REML.
- Added multiple grouped random-effect terms and correlated random coefficients.
- Added binomial-logit and Poisson-log GLMM fitting using PIRLS and a Laplace objective.
- Added prediction, simulation, covariance conversion, quadrature, PCA, and residual diagnostics.
- Bundled the modern Fortran `minqa` dependency for multidimensional BOBYQA optimization.
- Added numerical tests and independent `statsmodels` validation.
