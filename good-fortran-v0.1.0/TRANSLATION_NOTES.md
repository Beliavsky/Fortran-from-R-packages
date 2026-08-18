# Translation notes

## Source mapping

| R source | Fortran translation |
|---|---|
| `R/dgood.R` | `src/good_distribution.f90`: `dgood`, `good_logpmf` |
| `R/pgood.R` | `src/good_distribution.f90`: `pgood` |
| `R/qgood.R` | `src/good_distribution.f90`: `qgood` |
| `R/rgood.R` | `src/good_distribution.f90`: `rgood` |
| `R/good-internal.R::goodmean` | `src/good_distribution.f90`: `goodmean`, `good_moments` |
| `R/glm.good.R` | `src/good_glm.f90`: `glm_good` |
| `predict.glm.good` | `src/good_glm.f90`: `predict_good` |
| `summary.glm.good` | `src/good_glm.f90`: `summary_good` |

## Numerical changes

The R package obtains the normalizing polylogarithm through `copula::polylog` and switches to an asymptotic expression when the returned polylog overflows. The Fortran implementation instead accumulates the positive polylog series in scaled log-sum-exp form. This keeps the normalizer in log form and remains finite for cases such as `s=-170` where the unscaled polylog itself exceeds floating-point range.

The same series pass also computes expectations under the normalized weights. This gives stable analytic likelihood gradients and prediction gradients:

- `d log L / ds` uses `E(log N)`.
- `d log L / d eta` uses `N-E(N)` and the selected link derivative.
- The Good mean is `E(N)-1`.
- The derivative of the mean with respect to `s` is `-Cov(N,log N)`.
- The derivative with respect to the linear predictor uses `Var(N)`.

`maxLik`/`nlm` are replaced by a self-contained feasibility-preserving BFGS optimizer. The observed information matrix is obtained by finite-differencing the analytic negative-log-likelihood gradient and is inverted with pivoted Gauss-Jordan elimination.

The Fortran API accepts a design matrix directly instead of an R formula/data frame.

`rgood` uses inverse-CDF generation directly. The original `th` argument is retained for API familiarity and tail-search setup, but the generated variates are not deliberately truncated; this avoids the rare boundary artifacts in the R implementation.
