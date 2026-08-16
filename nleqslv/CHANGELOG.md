# Changelog

## 0.1.1

- Fixed strict FPM builds using `-Werror=implicit-interface`.
- Added explicit BLAS/LAPACK interfaces for all routines used by the numerical kernels.
- Collected the converted legacy kernels into `nleqslv_legacy_kernels`, giving all solver-to-solver calls explicit module interfaces.
- Added explicit callback interfaces for function and Jacobian procedures inside the legacy kernel layer.
- Revalidated all tests under `-Wimplicit-interface -Werror=implicit-interface` with runtime checking and at `-O2`.

## 0.1.0

- Initial standalone Fortran/FPM port of nleqslv 3.3.7.
- Converted upstream fixed-form Fortran numerical kernels to free-form F2018.
- Removed R and C runtime/interface requirements.
- Added modern typed callback API with option/result derived types.
- Added Newton and Broyden solver selection.
- Retained all seven upstream globalization modes.
- Added numerical, banded numerical, and user-supplied Jacobian support.
- Added automatic/fixed variable scaling and singular-Jacobian option.
- Added `search_zeros` and `test_nleqslv` computational helpers.
- Added FPM manifest, tests, example, and provenance documentation.
