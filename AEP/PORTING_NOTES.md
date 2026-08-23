# Porting notes

## Upstream coverage

Translated numerical routines from `R/AEP.R`:

| R | Fortran |
|---|---|
| `daep` | `daep` |
| `paep` | `paep` |
| `qaep` | `qaep` |
| `raep` | `raep` |
| `fitaep` | `fitaep` + `aep_fit_result` |
| `regaep` | `regaep` + `aep_reg_result` |

`welcome` and R-specific presentation/data handling are intentionally omitted.

## Random generation

Upstream `raep` uses a positive-stable mixture representation. The Fortran port generates exact AEP variates by inverse transformation through `qaep`. Both target the same distribution; inverse transformation removes the need for a separate stable-law rejection sampler.

## Estimation

The weighted fixed-point/EM-style updates from upstream are retained. R's `optimize`, `uniroot`, `solve`, `pgamma`, and `qgamma` are replaced by self-contained Fortran implementations.

## Source corrections

1. Upstream `fitaep` sets `n.p <- 3` despite estimating `alpha`, `sigma`, `mu`, and `epsilon`. The Fortran information criteria count four fitted parameters.
2. Upstream `regaep` uses `(S.T-S.E)/(p-1)*(n-p)*S.E` for its F statistic. The Fortran port uses the standard `((S.T-S.E)/(p-1))/(S.E/(n-p))`.
3. Upstream regression initialization effectively invokes `lm` on only the second column in a way that is not general for multiple predictors. The Fortran initialization uses ordinary least squares over all supplied predictors.
