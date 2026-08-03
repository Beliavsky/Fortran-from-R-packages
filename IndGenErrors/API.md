# API

Use the public module:

```fortran
use indgenerrors
```

All real inputs use `real(dp)`.

## `crosscor_2series(x, y, lag)`

Returns `type(lag_test_result)` with:

- `stat`: cross-correlations for lags `-lag:lag`
- `lags(1,:)`: lag labels
- `aggregate`: `n * sum(stat**2)`
- `p_aggregate`: chi-square upper-tail probability
- `n`, `status`

## `crosscor_3series(x, y, z, lag2, lag3)`

Returns `type(four_lag_test_result)` with components `xy`, `xz`, `yz`, and
`xyz`. The `xyz%lags` array has two rows containing the lags of the second
and third series. `aggregate` and `p_aggregate` combine all four subsets.

## `crossdep_2series(x, y, lag)`

Returns `type(dependence_two_result)` with `spearman`, `vdw`, and `savage`
components, each a `lag_test_result`.

## `crossdep_3series(x, y, z, lag2, lag3)`

Returns `type(dependence_three_result)`. Each of `spearman`, `vdw`, and
`savage` is a `four_lag_test_result`.

## `cvm_2series(x, y, lag)`

Returns `type(cvm_test_result)` with:

- `cvm`: normalized Cramer-von Mises statistics
- `p_cvm`: finite-sample empirical probabilities
- `wstat`, `p_wstat`: bias-corrected sum and Edgeworth probability
- `fstat`, `p_fstat`: Fisher combination and Edgeworth probability
- `lags`, `n`, `status`

## `cvm_3series(x, y, z, lag2, lag3)`

Returns `type(cvm_three_result)` with `xy`, `xz`, `yz`, and `xyz`
components, plus combined `wstat`, `fstat`, `p_wstat`, and `p_fstat`.

## Status codes

- `indgen_success = 0`
- `indgen_invalid_argument = 1`
- `indgen_numerical_error = 2`
