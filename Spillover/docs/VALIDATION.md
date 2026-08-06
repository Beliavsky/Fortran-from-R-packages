# Validation

The package contains five deterministic test programs.

## `test_var`

Simulates a two-variable VAR(1), fits it with a constant, checks coefficient and intercept recovery, verifies positive residual variances, and checks the first three MA coefficient matrices.

## `test_generalized_fevd`

Uses a zero-dynamics two-variable model with covariance

```text
[1.0  0.5]
[0.5  1.0]
```

and verifies the exact normalized generalized FEVD rows `(0.8, 0.2)`, the unnormalized cross share `0.25`, row normalization at every horizon, and total connectedness of 20 percent.

## `test_orthogonalized`

Checks the exact single-order Cholesky decomposition, the two-permutation average table with diagonal 87.5 and off-diagonal 12.5, total connectedness of 12.5 percent, compatibility-table layout, and the optional upstream `partial`/`total` label swap.

## `test_dynamic`

Simulates a two-regime VAR whose second regime introduces a strong cross-lag. The mean rolling connectedness increases from about 2.78 percent to about 12.08 percent. It also checks that net connectedness sums to zero in every window.

## `test_rolling_and_errors`

Checks rolling total and net dimensions, net balance, invalid-window handling, and the exact-permutation dimension guard.

Both `make check` and `make release` compile and run all tests. The checked build uses `-fcheck=all`; the release build uses `-O3`. Both use `-Wall -Wextra -Werror`.
