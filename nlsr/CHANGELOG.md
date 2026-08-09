# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of nlsr 2026.4.29 computational core.
- Added stabilized `nlfb` solver with bounds, masks, weights, damping, backtracking, and roff stopping.
- Added analytic and numerical residual-Jacobian interfaces.
- Added `jafwd`, `jaback`, `jacentral`, and Richardson `jand` equivalents.
- Added covariance/standard-error utilities.
- Added `SSlogisJN` numeric value/Jacobian/initializer utilities.
- Added six strict regression tests and two examples.
