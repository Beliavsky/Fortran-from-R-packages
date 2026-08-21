# Porting map

This file records how the computational pieces of `leaps` 3.2 map into the
Fortran package.

| Upstream component | Fortran translation |
| --- | --- |
| `src/leaps.f`: QR/subset algorithms | `src/leaps_lsq.f90`, `src/leaps_find_subsets.f90` using the supplied Alan Miller free-format modules |
| `src/leapshdr.f: MAKEQR` | `fit_qr` inside `src/leaps.f90`, using `startup` + repeated `includ` |
| `leaps.setup` | `regsubsets_fit` setup, weighting, column ordering, singularity handling, forced variables |
| `leaps.exhaustive` | `xhaust` through `regsubsets_fit(method='exhaustive')` |
| `leaps.backward` | `bakwrd` through `regsubsets_fit(method='backward')` |
| `leaps.forward` | `forwrd` through `regsubsets_fit(method='forward')` |
| `leaps.seqrep` | `seqrep` through `regsubsets_fit(method='seqrep')` |
| `summary.regsubsets` | result fields `rss`, `rsq`, `adjr2`, `cp`, `bic`, `model`, `valid` |
| `coef.regsubsets` | `model_coefficients` |
| `vcov.regsubsets` | optional `vcov` result of `model_coefficients` |
| `leaps` compatibility function | covered by exhaustive results plus the metric fields in `regsubsets_result` |
| formula/S3/printing code | omitted as R-specific infrastructure |
| `leaps.from.biglm` object adapter | omitted as R/`biglm`-specific infrastructure; the underlying QR/search algorithms are present |
| `plot.regsubsets` | omitted as requested |

## Compatibility notes

The metric formulas follow the R implementation in `R/leaps.R`, including
its no-intercept definition of `nullrss = sum(y**2)`.

For backward and forward selection with `nbest=1`, the R package defaults to
a nested path.  The Fortran wrapper preserves that behavior.  With
`nbest>1`, Miller's reporting machinery records multiple candidate subsets.

Predictor identifiers in the public result always refer to the original
columns of `x`, even when the setup phase internally reorders variables to
put forced-in variables first and forced-out/dependent variables last.
