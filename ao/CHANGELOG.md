# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `ao` 1.2.3.
- Ported alternating-optimization orchestration, partition generation,
  stopping logic, histories, best-process merging, and target splitting.
- Added standalone bounded BFGS, Nelder-Mead, and Newton block solvers for the
  external `optimizeR` base-optimizer boundary.
- Added strict GNU Fortran regression tests and examples.
