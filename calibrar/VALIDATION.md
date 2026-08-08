# Validation

Compiler used for release validation:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

Strict flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Release tests:

1. `test_gradient` -- forward/backward/central/Richardson gradients.
2. `test_fitness` -- calibration error/likelihood functions.
3. `test_random_kernel` -- truncated bounds, multivariate normal density, Gaussian grid.
4. `test_spline_objective` -- spline interpolation and weighted calibration terms.
5. `test_optimizers` -- BFGS, CG, Nelder-Mead, Hooke-Jeeves, SPG, and AHR-ES dispatch.
6. `test_calibrate` -- phased parameter activation and per-phase replicate counts.
7. `test_ahres` -- bounded multi-objective AHR-ES.
8. `test_stopping` -- calibrar stopping criteria.

Both examples are also compiled and run under the same flags.

The release archive is finally extracted into a new directory and rebuilt/tested from only the archived contents.
