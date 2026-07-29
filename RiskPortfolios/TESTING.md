# Testing

## Compiler and libraries

Validated with:

- GNU Fortran 14.2.0
- System BLAS
- System LAPACK

## Strict test flags

```text
-std=f2018
-Wall
-Wextra
-Werror
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-fbacktrace
-ffree-line-length-none
```

The optimized demonstration was also compiled with `-O2` and warnings treated
as errors.

## Test programs

`test/test_riskportfolios.f90` covers:

- All four mean estimators
- Both semideviation estimators
- All ten covariance estimators
- Matrix dimensions, symmetry, finite values, positive diagonals, and positive
  definiteness on a full-rank synthetic sample
- All seven portfolio methods
- Budget constraints
- Long-only constraints
- User lower and upper bounds
- Gross-exposure constraints
- Equal-risk-contribution accuracy

`test/test_reference_values.f90` compares the Fortran output against
independently computed fixed numerical references for:

- Naive, EWMA, and Bayes-Stein means
- Naive semideviation
- Sample covariance
- EWMA covariance
- Constant-correlation covariance
- Ledoit-Wolf market shrinkage
- Constant-correlation shrinkage
- Diagonal shrinkage
- One-parameter shrinkage
- Large-dimensional market shrinkage
- Bayes-Stein covariance

The factor estimator is tested structurally because the Fortran port uses the
principal-component factor replacement documented in `PORTING.md` rather than
R's `MASS::factanal` optimizer.

## Commands used

With FPM:

```text
fpm test
fpm run demo_riskportfolios
```

The validation environment did not contain the `fpm` executable. The manifest
was parsed independently, and all library, application, and test targets were
compiled directly with GNU Fortran using the strict flags above and
`-llapack -lblas`.
