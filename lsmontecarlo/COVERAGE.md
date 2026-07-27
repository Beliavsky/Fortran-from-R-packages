# Computational coverage

## Exported numerical routines

| Original R routine | Modern Fortran routine | Status |
|---|---|---|
| `AmerPutLSM` | `amer_put_lsm`, `american_put_lsmc` | Implemented |
| `AmerPutLSM_AV` | `amer_put_lsm_av`, `american_put_lsmc_antithetic` | Implemented |
| `AmerPutLSM_CV` | `amer_put_lsm_cv`, `american_put_lsmc_control` | Implemented |
| `AmerPutLSMPriceSurf` | `amer_put_lsm_price_surface` | Implemented |
| `AsianAmerPutLSM` | `asian_amer_put_lsm`, `asian_american_put_lsmc` | Implemented |
| `AsianAmerPutLSMPriceSurf` | `asian_amer_put_lsm_price_surface` | Implemented |
| `QuantoAmerPutLSM` | `quanto_amer_put_lsm`, `quanto_american_put_lsmc` | Implemented |
| `QuantoAmerPutLSM_AV` | `quanto_amer_put_lsm_av`, `quanto_american_put_lsmc_antithetic` | Implemented |
| `QuantoAmerPutLSMPriceSurf` | `quanto_amer_put_lsm_price_surface` | Implemented |
| `EuCallBS` | `eu_call_bs` | Implemented |
| `EuPutBS` | `eu_put_bs` | Implemented |
| `fastGBM` | `fast_gbm` | Implemented |
| `firstValueRow` | `first_value_row` | Implemented |
| `price` | generic `price` | Implemented |

Price-surface functions return a typed `price_surface` object containing the
volatility vector, strike vector, and numerical price matrix.

## R object and presentation routines

The following R-specific facilities are not compiled:

- `plot.PriceSurface`, which calls R's `persp`
- S3 `print.*` methods
- S3 `summary.*` methods

Their numerical data are available directly through `option_result` and
`price_surface`. Surface summary helpers `surface_mean`, `surface_minimum`, and
`surface_maximum` are included.

## Internal replacements

The original dependencies are replaced as follows:

| R dependency | Fortran replacement |
|---|---|
| `rnorm` | native Box-Muller normal generator |
| `rmvnorm` | native correlated-normal generator |
| `lm` | scaled rank-revealing QR least squares |
| `pnorm` | intrinsic `erfc` normal CDF |
| matrix/list/S3 infrastructure | allocatable arrays and typed derived types |

No external numerical or statistical library is required.
