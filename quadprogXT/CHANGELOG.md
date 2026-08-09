# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `quadprogXT` 0.0.6.
- Added `build_qp_xt`, `solve_qp_xt`, `normalize_constraints`, and
  `convert_to_compact`.
- Vendored the supplied `quadprog-fortran` translation as an FPM dependency.
- Added absolute-value and absolute-change constraints and objective penalties.
- Added dense/compact, normalization, and factorized-QP paths.
- Added strict GNU Fortran regression tests and examples.
