# Computational coverage

## Translated exports

| Upstream routine | Fortran interface | Status |
|---|---|---|
| `dvasicek` | `vasicek_density` | Complete |
| `pvasicek` | `vasicek_cdf` | Complete |
| `qvasicek` | `vasicek_quantile` | Complete |
| `rvasicek` | `random_vasicek` | Complete |
| `.vasicek_mu` | `effective_probit_mean` | Complete |
| `vasicekfit` | `fit_vasicek` | Complete numerical translation |
| `coef.vasicekfit` | `coefficients` | Complete |
| `fitted.vasicekfit` | `fit%fitted` | Complete |
| `residuals.vasicekfit` | `fit%residuals` | Complete |
| `predict.vasicekfit`, link | `predict_link` | Complete |
| `predict.vasicekfit`, mean response | `predict_response` | Complete |
| `predict.vasicekfit`, quantiles | `predict_quantiles` | Complete |
| `vcov.vasicekfit`, IID | `vasicek_covariance` | Complete |
| `vcov.vasicekfit`, HAC | `vasicek_covariance(use_hac=.true.)` | Methodologically complete; native Bartlett HAC |
| `confint.vasicekfit` | `vasicek_confidence_intervals` | Complete |

## Internal numerical support

- Normal density, CDF, and inverse CDF
- Seedable standard-normal random generation
- OLS regression and covariance
- Matrix inversion with pivoting
- Sample variance
- Bartlett/Newey-West long-run covariance
- Exact matrix symmetrization

## R-specific exclusions

The following components have no compiled numerical counterpart:

- R formulas, model frames, terms objects, and variable names
- S3 print and summary formatting
- `printCoefmat` output
- `sandwich` package configuration and prewhitening adapters
- R documentation and testthat infrastructure

Fortran callers provide the response vector and predictor matrix directly.
