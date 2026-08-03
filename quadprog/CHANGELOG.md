# Changelog

## 1.5.8-fortran.1

- Modernized both dense and compact Goldfarb-Idnani solvers to free-form
  Fortran 2018 modules.
- Added typed FPM API and `qp_result` diagnostics.
- Replaced external LINPACK/BLAS dependencies with self-contained routines.
- Added strict regression, factorized-input, compact-storage, equality,
  infeasibility, positive-definiteness, and KKT tests.
- Added examples, demo, documentation, and cross-platform validation scripts.
