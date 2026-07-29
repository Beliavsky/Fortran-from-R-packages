# Testing

## Automated coverage

`test/test_risk.f90` checks:

- density, CDF, and quantile values for all six built-in distributions
- Student t CDF/quantile inversion
- a user-defined triangular callback distribution
- all 26 translated risk-measure procedures
- scalar and vector generic overloads
- finite and infinite integration intervals
- source-compatible expectile behavior
- invalid-domain NaN handling

Reference cases include analytic normal and uniform results. An additional
independent SciPy comparison used a shifted and scaled normal distribution and
a nonstandard uniform distribution. Across eight representative calculations,
the largest absolute Fortran-versus-SciPy difference was approximately
`6.3e-12`.

## Compiler and flags

The final test run used GNU Fortran 14.2.0 with:

```text
-std=f2018
-Wall
-Wextra
-Werror
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-fbacktrace
-O2
```

The result was:

```text
All Risk tests passed.
```

## Commands

With FPM:

```text
fpm test
fpm run
fpm run --example custom_distribution
```

Direct GNU Fortran compilation:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -fbacktrace -O2 \
  src/risk_kinds.f90 src/risk_math.f90 \
  src/risk_distributions.f90 src/risk_numerics.f90 \
  src/risk_measures.f90 src/risk.f90 test/test_risk.f90 \
  -o test_risk

./test_risk
```
