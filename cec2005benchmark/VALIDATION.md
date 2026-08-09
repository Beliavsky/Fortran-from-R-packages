# Validation

Compiler used during translation validation:

```text
GNU Fortran 14.2.0
```

Strict flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

## Regression tests

### `test_reference_data`

For each of the 25 functions, reads the ten 50-dimensional vectors and ten
expected values from the upstream `test_data_func*.txt` files.  Noise is
disabled, matching upstream `tests/test_data.R`.

Total comparisons: **250**.

Worst relative error observed: approximately **2.35e-14**.

### `test_global_optima`

Reproduces upstream `tests/test_global_optima.R` for every function and all
four dimensions (2, 10, 30, 50), including the F5, F8, and F20 boundary
modifications.

Total optimum checks: **100**.

Worst relative error observed: approximately **2.35e-14**.

### `test_context_batch_noise`

Checks persistent-context versus convenience batch evaluation and verifies the
nonnegative multiplicative noise behavior of F4.

### `test_validation`

Checks invalid function ids, unsupported dimensions, and input-size errors.

## Examples

- `basic`
- `reuse_context`

Both compile and run under the same strict flags.
