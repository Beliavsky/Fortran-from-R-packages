# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64
- No external numerical libraries
- FPM was not installed in the validation environment

The release contains 1780 lines across 17 Fortran source,
test, demo, and example files.

## Compiler configurations

Checked build:


- -std=f2018
- -O0 -g
- -Wall -Wextra -Wconversion-extra -Werror
- -fcheck=all -fbacktrace
- -ffree-line-length-132

Optimized build:


- -std=f2018
- -O2
- -Wall -Wextra -Wconversion-extra -Werror
- -ffree-line-length-132

Both configurations passed all tests, the demo, and both examples.

## Test results


test_conversions: PASS
test_models: PASS
test_pca: PASS
test_risk_analysis: PASS
validation (debug): PASS
validation (optimized): PASS

## Numerical references

The validation suite includes:

- Exact synthetic Nelson-Siegel parameter recovery.
- Exact synthetic Svensson parameter recovery.
- Natural and FMM/not-a-knot spline values independently calculated with
  SciPy's CubicSpline implementation.
- Fixed par-to-zero bootstrap values and annual/semi-annual round trips.
- Hand-calculated coupon-bond price, duration, and convexity values.
- Fixed Z-spread and key-rate-duration identities.
- Flat-curve carry and roll-down identities.
- Independent NumPy covariance-eigenvalue PCA variance shares.
- Orthonormal PCA loading checks and scaled-PCA smoke tests.

## Release audits

The release audit checks:

- Valid TOML syntax for fpm.toml.
- ASCII-only translated source and documentation.
- Free-form Fortran lines no longer than 132 columns.
- SPDX MIT headers in every Fortran file.
- implicit none in every Fortran program and module.
- No object files, module files, executables, or build directories in the ZIP.
- SHA-256 manifests for the supplied archive, unmodified original source, and
  translated release files.

The exact final ZIP is extracted into a clean directory and both validation
configurations are rerun before delivery.
