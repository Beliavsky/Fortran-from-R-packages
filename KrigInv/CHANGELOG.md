# Changelog

## 0.2.0

- Added a vendored GPL-3-compatible subset of DiceKriging-fortran 0.1.0.
- Added `fit_krig_model` for MLE/PMLE/LOO Gaussian-process fitting.
- Added `init_dice_krig_model` for fixed parameters with DiceKriging covariance conventions.
- Added genuine covariance/variance re-estimation to `update_krig_model`.
- Reproduced KrigInv's default `CovReEstimate = model@param.estim` behavior.
- Added `cov_reestimate`, `trend_reestimate`, and `nugget_reestimate` controls to sequential EGI updates.
- Added DiceKriging-backed SK/UK prediction and posterior covariance.
- Added standard DiceKriging interaction and quadratic trend bases.
- Added power-exponential, isotropic, and nonlinear scaling covariance support through the fitted-model backend.
- Added deterministic/noisy re-estimation and EGI regression tests.
- Preserved the v0.1.0 fixed-parameter `init_krig_model` path for backward compatibility.

## 0.1.0

- Initial modern-Fortran/FPM translation of KrigInv's non-plotting computational core.
