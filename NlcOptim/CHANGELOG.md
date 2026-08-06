# Changelog

## 0.1.1

- Fixed `fpm build` failure caused by duplicate `demo_problem` module and `demo_nlcoptim` program names in `app/` and `example/`.
- Renamed the example units to `basic_nonlinear_problem` and `basic_nonlinear`.

## 0.1.0

- Initial modern Fortran translation of NlcOptim 0.6.
- Added typed SQP options, results, multipliers, and status codes.
- Added finite-difference and analytic derivative interfaces.
- Integrated the supplied modern Fortran `quadprog` dependency.
- Added damped BFGS, exact-penalty line search, and elastic feasibility QP.
- Added FPM/Make builds, demonstration, five tests, and provenance documents.
