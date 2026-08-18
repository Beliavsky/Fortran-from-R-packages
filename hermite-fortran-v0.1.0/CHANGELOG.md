# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `hermite` 1.2.1.
- Ported exact and approximate generalized-Hermite distribution routines.
- Added exact log-sum PMF for validation and approximation fallback.
- Ported the raw numerical `glm.hermite` likelihood and fitting workflow.
- Added fixed/automatic order fitting, covariance/Hessian, LR test, and AIC.
- Replaced the upstream unrelated standard-normal PMF fallback with exact mass.
