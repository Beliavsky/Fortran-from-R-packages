# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs:

1. `test_qp_example` — upstream simplex/nonnegative QP example;
2. `test_isotonic` — isotonic projection;
3. `test_factorized` — factorized diagonal input and source-compatible output;
4. `test_semidefinite` — rank-deficient PSD quadratic matrix;
5. `test_status` — iteration-limit/nonconvergence bookkeeping.

Examples:

- `simplex_qp`;
- `isotonic`.

The release archive is also extracted into a fresh directory and rebuilt from
only its contents before publication.
