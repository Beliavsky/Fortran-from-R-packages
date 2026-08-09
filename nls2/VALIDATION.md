# Validation

Validation compiler:

```text
GNU Fortran 14.2.0
```

Strict flags:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Regression programs cover:

1. nonlinear exponential fitting with numerical and analytical Jacobians and weights;
2. source-compatible brute-force grid sizing/order and default multi-start fitting;
3. Latin-hypercube stratification and bounds;
4. partially-linear variable-projection fitting and brute-force evaluation;
5. singular-Jacobian start evaluation and bounded `port` compatibility.

Both examples are also compiled/run by `scripts/test_gfortran.sh` and `scripts/test_gfortran.bat`.

The final release is additionally validated after extraction into a new directory so the build cannot depend on untracked working files.
