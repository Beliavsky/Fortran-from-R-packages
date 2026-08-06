# Changelog

## 0.1.1

- Rename each test helper module uniquely so FPM does not reject the package for duplicate `test_support` module names.
- No numerical library code or public API changed.

## 0.1.0

- Translated all exported computational StReg models to modern Fortran.
- Added exact block-Toeplitz covariance construction.
- Added self-contained linear algebra, BFGS/Nelder-Mead optimization,
  Student-t probabilities, Anderson-Darling diagnostics, and delta-method
  inference.
- Added FPM and Make builds, examples, tests, API mapping, provenance, and
  source-compatibility documentation.
