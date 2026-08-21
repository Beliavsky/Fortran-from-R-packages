# Validation

Validation was performed with GNU Fortran 14.2.0.

Compiler options:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -fcheck=all -fbacktrace
```

The validation build starts from an empty object/module directory and compiles the same `src/`, `test/`, and `example/` units referenced by the FPM project.

## Tests

```text
test_core:    PASS
test_fit:     PASS
test_rng_gof: PASS
```

`test_core` checks fixed PMF reference values, total probability, reduction to independent Poisson variables at `lambda3=0`, probability-grid bounds, and contingency-table aggregation.

`test_fit` uses a fixed 80-observation sample and independent SciPy 1.17 calculations. Reference targets include:

```text
lambda3 MLE       2.48070493735253
lambda1 MLE       2.38179506264747
lambda2 MLE       4.31929506264747
fitted logLik    -355.3619376814744
independent logLik -363.6849481214147
profile-grid CI   [1.39, 3.30]
```

The full fit also checks the observed profile-curvature variance and the translated Kawamura asymptotic variance against independently calculated values.

`test_rng_gof` checks simulated marginal means/covariance for both small and large Poisson rates, exercises the PTRS large-rate branch, and verifies the Monte Carlo GOF APIs.

The demo produces, for its seeded sample:

```text
true lambda:          3.000000    5.000000    2.000000
estimated lambda:     2.571448    4.718115    2.258552
P(X=3,Y=4):         0.01754045
```

FPM itself was not installed in the validation environment. `fpm.toml` was parsed with Python's TOML parser and the exact FPM source/test units were compiled manually with `gfortran`.
