# API map

| Upstream R routine | Fortran counterpart | Notes |
|---|---|---|
| `dskellam` | generic `dskellam` | Scalar or integer-vector input; scalar rates |
| `pskellam` | generic `pskellam` | Lower/upper tail and log-probability options |
| `qskellam` | generic `qskellam` | Returns `integer(i8)` quantiles |
| `rskellam` | `rskellam` | Returns an allocatable `integer(i8)` vector |
| `dskellam.sp` | `dskellam_sp` | Correctly honors the log option |
| `pskellam.sp` | `pskellam_sp` | Lugannani-Rice approximation with continuity correction |
| `skellam.mle` | `fit_skellam_mle` | Same one-dimensional mean constraint; transformed optimization |
| `skellam.reg` | `fit_skellam_regression` | Matrix-first API; optional intercept and response ordering |

Additional Fortran routines include `skellam_log_pmf`,
`skellam_log_likelihood`, moment functions, and `seed_random_number`.

## Result types

`skellam_mle_result` contains rates, standard errors, covariance, likelihood,
iteration counts, status, and convergence state.

`skellam_regression_result` contains both coefficient vectors, covariance,
standard errors, Wald statistics, p-values, fitted rates, likelihood, and
optimization diagnostics.
