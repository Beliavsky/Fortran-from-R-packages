# Build report

- Build date: 2026-08-03
- Compiler: GNU Fortran (Debian 14.2.0-19) 14.2.0
- Input archive SHA-256: `86b5b597c3b0a6b76011b2989ba36323240953e49ff38397784a086c49b899f9`
- FPM manifest: parsed successfully with Python's TOML parser
- FPM executable: not installed in the execution environment; the equivalent
  library, test, example, and application sources were compiled directly with
  GNU Fortran.

## Checked build

Command family: GNU Fortran 2018 with warnings as errors, implicit typing and
implicit externals disabled, runtime bounds/type checks, and backtraces.

Results:

- `test_statistics`: PASS
- `test_covariance_estimators`: PASS
- `test_pca`: PASS
- `test_discriminant`: PASS
- `test_multivariate_tests`: PASS
- `demo_rrcov`: PASS
- `example_covariance`: PASS

## Optimized build

GNU Fortran `-O3` build: all five tests PASS. Two optimizer-only false-positive
warnings concerning unallocated derived-type descriptors are explicitly
suppressed in `scripts/build_optimized.sh`; all other warnings remain errors.

## FPM-like default build

The source and all tests also compile and run with ordinary gfortran defaults,
without the GNU free-line-length extension or special module-order workarounds.
All Fortran source lines are at most 132 characters.
