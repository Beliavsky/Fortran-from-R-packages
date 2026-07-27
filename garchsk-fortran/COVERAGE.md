# Computational coverage

Source package: `GARCHSK` 0.1.0.

| Upstream routine | Fortran routine | Status |
|---|---|---|
| `skewness` | `skewness` | Complete |
| `kurtosis` | `kurtosis` | Complete |
| `garchsk_construct` | `garchsk_construct` | Complete |
| `garchsk_lik` | `garchsk_lik` | Complete |
| `garchsk_ineqfun` | `garchsk_ineqfun` | Complete |
| `garchsk_est` | `garchsk_est` | Complete; native optimizer |
| `garchsk_fcst` | `garchsk_fcst` | Complete; corrected recursion |
| `gjrsk_construct` | `gjrsk_construct` | Complete |
| `gjrsk_lik` | `gjrsk_lik` | Complete |
| `gjrsk_ineqfun` | `gjrsk_ineqfun` | Complete |
| `gjrsk_est` | `gjrsk_est` | Complete; native optimizer |
| `gjrsk_fcst` | `gjrsk_fcst` | Complete; corrected recursion |

## Supporting numerical algorithms

The translation additionally provides:

- sample mean, variance, covariance, skewness, and Pearson kurtosis
- constrained Nelder-Mead optimization
- numerical Hessians
- pivoted Gauss-Jordan matrix inversion
- corrected inverse-Hessian covariance and standard errors
- typed paths, forecasts, and estimation results

## Not compiled

- The bundled `GBP.rda` Bloomberg data object
- R package documentation and namespace infrastructure
- `Rsolnp` bindings

The original data and source remain under `original/GARCHSK-0.1.0` for
provenance. No numerical function was excluded.
