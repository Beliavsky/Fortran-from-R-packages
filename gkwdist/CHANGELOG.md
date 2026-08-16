# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of gkwdist 1.1.4.
- Ported all 50 non-pipe computational exports.
- Added standalone GKw/BKw/KKw/EKw/McDonald/Kumaraswamy/Beta d/p/q/r routines.
- Added negative log-likelihoods and analytical gradients/Hessians.
- Added method-of-moments starting-value estimation with Nelder-Mead.
- Replaced R Mathlib beta/special-function calls with standalone Fortran.
- Added tests against closed-form references, nesting identities, numerical
  differentiation, RNG properties, and API coverage.
