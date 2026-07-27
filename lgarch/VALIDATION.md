# Validation

## Toolchain

- GNU Fortran 14.2.0
- LAPACK and BLAS from the system libraries
- Fortran 2018 mode

## Runtime-checked build

Flags:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -fcheck=all -fbacktrace
```

The following passed:

- GPL-2.0-only source-header and metadata audit.
- Vector and matrix lag/difference reference tests.
- Deterministic arbitrary-order univariate simulation recursion.
- Univariate ARMA recursion with zero-return conditional imputation.
- Direct Gaussian ARMA likelihood reference calculation.
- LS, Gaussian QML, and CEX2 estimation paths.
- Mean-correction path.
- Exogenous-regressor path.
- Custom initial values, bounds, and objective penalties.
- Numerical Hessian and covariance allocation.
- Multivariate Normal mean and covariance simulation checks.
- Deterministic multivariate log-GARCH simulation recursion.
- Multivariate zero-return recursion.
- Direct multivariate Gaussian likelihood reference calculation.
- Full CCC log-GARCH(1,1) fit path.
- Multivariate custom initial values, bounds, and penalty path.
- Demonstration, exogenous-regressor example, and both CSV modes.

## Optimized build

Flags:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -fbacktrace
```

The same tests and programs passed.

## Interpretation

The deterministic tests validate recursion orientation, parameter packing,
zero-return behavior, and likelihood formulas. Simulation moment tests use
finite tolerances. Estimation tests verify executable end-to-end workflows,
finite valid results, positive fitted standard deviations, parameter and
covariance dimensions, and custom-control propagation.

Because the translation uses a different optimizer and RNG, exact R random
streams and exact optimizer endpoints are not claimed. R was not available in
the validation environment, so direct run-by-run R comparison is not claimed.
