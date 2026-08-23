# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of RSpectra 0.16-2 computational interfaces.
- Added ARPACK-backed symmetric and general partial eigensolvers.
- Added dense real symmetric shift-and-invert.
- Added exact dense transformed-spectrum selection for shifted general matrices, including complex shifts.
- Added native CSR sparse operator support.
- Added extensible user-defined `linear_operator` API.
- Added truncated SVD for tall/wide matrices with center/scale support.
- Added LAPACK full-spectrum/full-SVD fallbacks.
- Added strict numerical test suite and examples.
