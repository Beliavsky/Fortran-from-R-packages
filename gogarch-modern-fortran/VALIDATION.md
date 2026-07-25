# Validation

Validation was performed on 2026-07-23 with GNU Fortran 14.2.0, BLAS, and
LAPACK.

## Runtime-checked build

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace
```

## Optimized build

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
-Werror -fbacktrace
```

Both builds passed:

```text
GPL-2.0-or-later source license checks passed.
Core numerical tests passed.
Higher-order GARCH, APARCH, and distribution tests passed.
GO-GARCH estimator and workflow tests passed.
Extended GO-GARCH specification tests passed.
```

The demo executed successfully. The CSV application executed successfully for
all four estimators with the default model, for Student GARCH(2,1), and for
skew-Student APARCH(1,1,1).

## What the tests establish

- all six conditional densities numerically integrate to approximately one;
- random generators have sample mean approximately zero and variance
  approximately one under fixed seeds;
- every distribution is used in a likelihood fit;
- GARCH(2,2) and APARCH recursions, fitting, and forecasts execute with finite
  positive conditional variances;
- delta, leverage, shape, and skew fitting paths execute;
- ICA, MM, NLS, and ML accept a Student GARCH(2,1) factor specification;
- an APARCH skew-Student GO-GARCH model fits, forecasts, and exposes complete
  parameter vectors;
- existing release-0.1.0 GARCH(1,1) and estimator tests continue to pass.

## Equivalence limits

The tests establish internal numerical consistency and executable coverage.
They do not establish parameter-for-parameter equivalence with R `fGarch`,
`fastICA`, `optim`, or `nlminb`. R and `fpm` were unavailable in the validation
environment, so neither R comparison tests nor an `fpm` build are claimed.
