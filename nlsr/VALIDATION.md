# Validation

Compiler used during translation validation:

```text
GNU Fortran 14.2.0
```

Strict flags:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Regression programs:

1. `test_hobbs` - bounded Hobbs weed nonlinear regression with analytic Jacobian.
2. `test_jacobians` - forward/backward/central/Richardson numerical Jacobians and masks.
3. `test_linear_weights` - weighted one-parameter nonlinear-LS path and covariance calculation.
4. `test_mask_bounds` - exact parameter masking using `lower == upper`.
5. `test_logistic` - `SSlogisJN`-style initialization plus nonlinear fit.
6. `test_all_masked` - source-compatible immediate return when every parameter is fixed.

Representative Hobbs result:

```text
b1 = 196.1862600154
b2 = 49.0916392072
b3 = 0.3135697308
RSS = 2.5872773953
```

Both examples are also compiled and run by the strict validation scripts.
