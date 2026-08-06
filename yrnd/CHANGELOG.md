# Changelog

## 0.1.0

- Initial modern Fortran translation of the yrnd 0.1.5 computational code.
- Added lognormal-mixture option calibration for European, American, and
  futures-style-margin options.
- Added risk-neutral future-price and STIR-rate distributions.
- Added bond cash-flow schedules, yield transformations, net-basis calculations,
  and CTD probabilities using the translated tvm dependency.
- Added Gaussian-copula yield-spread simulation and KDE output.
- Added FPM packaging, tests, examples, validation scripts, and source provenance.
