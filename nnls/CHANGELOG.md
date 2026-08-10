# Changelog

## 0.1.0

- Initial modern Fortran translation of `nnls` 1.6.
- Added NNLS and mixed-sign NNNPLS solvers.
- Replaced fixed-form/GOTO work-array implementation with module-based
  Lawson-Hanson active-set logic and Householder QR passive solves.
- Added FPM package, examples, strict regression tests, provenance, and license
  files.
