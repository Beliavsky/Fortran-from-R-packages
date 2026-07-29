# Changelog

## 0.1.1-fortran.1

- Translated the complete RM2006 covariance algorithm to modern Fortran.
- Added an FPM package manifest and standard source layout.
- Added `rm2006` and `rm2006_covariance` public APIs.
- Added public multiscale time-scale and weight calculation.
- Added status-code validation for dimensions, parameters, and finite data.
- Reduced temporary storage by computing outer products on demand.
- Corrected the original small-sample endpoint indexing edge case.
- Added deterministic numerical tests and two runnable examples.
- Preserved GPL-2.0-or-later licensing and original attribution.
