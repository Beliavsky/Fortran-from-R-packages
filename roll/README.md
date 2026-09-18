# roll

Modern free-form Fortran translation of the computational API of the R package [`roll`](https://github.com/jasonjfoster/roll), version 1.2.1.

The package provides deterministic rolling/expanding statistics for vectors and matrices, including weighted moments, type-2 quantiles, pairwise covariance/correlation/crossproducts, and rolling weighted linear regression.

## Build

Place this directory beside the shared packages `rfortran-core` and `rfortran-linalg` in the root of `Fortran-from-R-packages`, then run:

```text
fpm build
fpm test
fpm run --example basic_roll
```

`rfortran-core` supplies the shared `dp` real kind. `rfortran-linalg` supplies the checked linear solves used by `roll_lm`; no BLAS/LAPACK source or system-library link is copied into this package.

## Numerical API

Import the public `roll` module:

```fortran
use roll, only : dp, roll_mean, roll_cov, roll_lm, roll_lm_result
```

Most univariate procedures are generic over rank-1 and rank-2 inputs. Matrix observations are rows and variables are columns, matching the R package. Real missing values are IEEE quiet NaNs. Logical data use integer zero/one values and the public sentinel `roll_na_logical` for missing values.

For the pairwise functions, vector inputs return a vector and matrix inputs return a cube with shape `(ncol(x), ncol(y), nrow(x))`. When `y` is omitted, `x` is paired with itself. `roll_lm` returns a `roll_lm_result` containing `coefficients`, `r_squared`, and `std_error` arrays.

The Fortran `online` argument is accepted for compatibility but both values use the same deterministic direct-window implementation. Therefore the numerical target is preserved without reproducing RcppParallel scheduling or every online worker's floating-point operation order.

## Translation coverage

Package status: **substantial**. **18 of 18 (100%)** exported computational R functions are mapped. This fraction measures mapped computational functions, not complete R compatibility. R/xts/zoo attributes, names/dimnames, `.Call` wrappers, warning text, and threaded execution are intentionally omitted.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `roll_all` | `roll_univariate`, `roll_matrix` | `roll_all_vec`, `roll_all_mat` | substantial |
| `roll_any` | `roll_univariate`, `roll_matrix` | `roll_any_vec`, `roll_any_mat` | substantial |
| `roll_sum` | `roll_univariate`, `roll_matrix` | `roll_sum_vec`, `roll_sum_mat` | substantial |
| `roll_prod` | `roll_univariate`, `roll_matrix` | `roll_prod_vec`, `roll_prod_mat` | substantial |
| `roll_mean` | `roll_univariate`, `roll_matrix` | `roll_mean_vec`, `roll_mean_mat` | substantial |
| `roll_min` | `roll_univariate`, `roll_matrix` | `roll_min_vec`, `roll_min_mat` | substantial |
| `roll_max` | `roll_univariate`, `roll_matrix` | `roll_max_vec`, `roll_max_mat` | substantial |
| `roll_idxmin` | `roll_univariate`, `roll_matrix` | `roll_idxmin_vec`, `roll_idxmin_mat` | substantial |
| `roll_idxmax` | `roll_univariate`, `roll_matrix` | `roll_idxmax_vec`, `roll_idxmax_mat` | substantial |
| `roll_median` | `roll_univariate`, `roll_matrix` | `roll_median_vec`, `roll_median_mat` | substantial |
| `roll_quantile` | `roll_univariate`, `roll_matrix` | `roll_quantile_vec`, `roll_quantile_mat` | substantial |
| `roll_var` | `roll_univariate`, `roll_matrix` | `roll_var_vec`, `roll_var_mat` | substantial |
| `roll_sd` | `roll_univariate`, `roll_matrix` | `roll_sd_vec`, `roll_sd_mat` | substantial |
| `roll_scale` | `roll_univariate`, `roll_matrix` | `roll_scale_vec`, `roll_scale_mat` | substantial |
| `roll_cov` | `roll_multivariate` | `roll_cov_vec`, `roll_cov_mat` | substantial |
| `roll_cor` | `roll_multivariate` | `roll_cor_vec`, `roll_cor_mat` | substantial |
| `roll_crossprod` | `roll_multivariate` | `roll_crossprod_vec`, `roll_crossprod_mat` | substantial |
| `roll_lm` | `roll_multivariate` | `roll_lm_vec`, `roll_lm_mat` | substantial |

`fpm.toml` contains the same machine-readable coverage count and one mapping entry per function.

## Compatibility notes

- Weighted sums/products/means use the final weight for the current observation and progressively earlier weights for earlier lags, matching upstream alignment.
- Variance and covariance use the upstream unbiased weighted denominator `sum(w) - sum(w**2)/sum(w)`.
- Quantiles implement the inverse empirical distribution with averaging at discontinuities (Hyndman-Fan type 2), including the upstream weighted extension.
- `complete_obs=.true.` on matrix statistics drops a row from a window when any matrix column is missing; pairwise covariance/correlation/crossproducts use pairwise rows when false.
- `roll_lm` follows upstream complete-row regression semantics even when `complete_obs=.false.`; the R warning itself is omitted.
- The R package's RcppParallel online/offline algorithms are not duplicated. Direct-window evaluation supports arbitrary supplied weight vectors, including cases for which upstream would warn or reject `online=.true.`.
- R classes, xts/zoo indexes, dimnames, and list presentation are not represented in the numerical Fortran API.

## Validation

See `VALIDATION.md` for the exact commands and environment used for this translation.

## License and provenance

The upstream package is GPL (>= 2). This translation is distributed under GPL-2.0-or-later. See `LICENSE`, `NOTICE.md`, and `PROVENANCE.md`.
