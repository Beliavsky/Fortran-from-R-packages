# Validation

Validation environment:

- GNU Fortran 14.2.0
- System LAPACK and BLAS
- Linux x86-64

## Debug configuration

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -O0 -g -fcheck=all -fbacktrace
```

Command:

```sh
make debug-check
```

Results:

```text
GPL-2.0-or-later source license checks passed.
Univariate and nonlinear tests passed.
Multivariate linear and cointegration tests passed.
Threshold multivariate and nonlinear impulse-response tests passed.
debug build, tests, and applications passed.
```

## Optimized configuration

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -O2 -fbacktrace
```

Command:

```sh
make release-check
```

Results:

```text
GPL-2.0-or-later source license checks passed.
Univariate and nonlinear tests passed.
Multivariate linear and cointegration tests passed.
Threshold multivariate and nonlinear impulse-response tests passed.
release build, tests, and applications passed.
```

## Tested numerical coverage

The suites exercise simulation, fitting, forecasting, order/rank selection, threshold regime counts, local-search equivalence, residual/wild/block resampling, rolling forecasts, VAR/VECM representations, Johansen statistics, IRFs, FEVD, GIRFs, and nonlinear/unit-root statistics. They use fixed random seeds and check parameter recovery, dimensions, finiteness, probability ranges, regime occupancy, FEVD normalization, and direct-versus-box local-linear equivalence.

No claim of exact R optimizer endpoints, random streams, external-package results, or bootstrap critical values is made.
