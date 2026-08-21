# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of BGFD 0.1.
- Added all 16 Bell-G/complementary Bell-G distribution families.
- Added d/p/q/r/s/h operations and generic run-time family API.
- Added all 16 MLE/goodness-of-fit wrappers using the translated AdequacyModel
  numerical layer.
- Corrected the four upstream Weibull/exponentiated-Weibull inverse-CDF powers.
- Defined consistent log-survival and log-hazard semantics.
- Added strict numerical reference, identity, derivative, RNG, option, and fit
  tests.
