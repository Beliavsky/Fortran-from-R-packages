# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `lbfgsb3c`.
- Converted the L-BFGS-B 3.0 fixed-form source to a free-form module.
- Added native analytic-gradient and finite-difference callback APIs.
- Added typed controls/results, scalar/vector bounds, user data, monitoring,
  cancellation, non-finite checks, and `lbfgsb3c` parameter tolerances.
- Replaced the R-provided LINPACK dependency with self-contained `dpofa` and
  `dtrsl` equivalents.
- Added tests based on the package's bounds, Chebyquad, and generalized
  Rosenbrock examples.
