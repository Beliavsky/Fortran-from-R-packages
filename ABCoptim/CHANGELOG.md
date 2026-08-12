# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of ABCoptim 0.15.0.
- Added separate `abc_optim` and `abc_cpp` entry points reflecting the
  package's pure-R and Rcpp algorithm semantics.
- Added continuous box-constrained and R-style binary optimization.
- Added `parscale`, `fnscale`, employed/onlooker/scout phases, persistence
  stopping, fitness transformation, histories, and diagnostics.
- Added deterministic standalone RNG, tests, examples, and license/provenance
  material.
