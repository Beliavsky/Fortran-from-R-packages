# Validation

## Environment

- GNU Fortran 14.2.0
- LAPACK and BLAS from the Debian runtime
- Fortran 2018 mode

## Debug configuration

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace
```

## Optimized configuration

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
-Werror -fbacktrace
```

## Executed suites

1. Chaotic maps and RK4
   - Exact Tent, Logistic, Henon, and Ikeda recursions
   - One-step Lorenz and Rossler RK4 reference values
2. Embedding and nonlinear statistics
   - Regular/explicit/matrix embeddings
   - Exact correlation integral
   - Correlation-dimension monotonicity
   - Finite-sample mutual-information convention
   - Recurrence and space-time separation outputs
3. Neighbor and Lyapunov calculations
   - False-neighbor fractions
   - Direct/box neighbor index and distance equivalence
   - Reduced exact-distance evaluations for box search
   - Direct/box Lyapunov-path equivalence on Lorenz data
4. Hypothesis tests
   - Hand-derived BDS correlations, Dechert `k`, and normalized statistics
   - Deterministic White-test random-weight and regression reference values
   - Terasvirta polynomial expansion and strong nonlinear-series detection
   - Exact runs-test example
   - Generic test dispatcher

The demo, CSV analyzer, and nonlinear-test example were also compiled and run.
Every Fortran file was checked for its GPL-2.0-or-later SPDX identifier and
license notice.

`fpm` was not installed and is not claimed as tested.
