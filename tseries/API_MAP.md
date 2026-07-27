# API map

Original package: R `tseries` 0.10-62.

## Computational functions

| R function | Fortran procedure | Status |
|---|---|---|
| `adf.test` | `adf_test` | Translated; experimental |
| `arma` | `arma_fit` | Translated with internal Nelder-Mead CSS optimizer |
| `bds.test` | `bds_test` | Translated with direct low-memory pair counting |
| `garch` | `garch_fit` | Translated with constrained parameter transform |
| `jarque.bera.test` | `jarque_bera_test` | Translated |
| `kpss.test` | `kpss_test` | Translated |
| `maxdrawdown` | `maximum_drawdown` | Translated; returns the first maximum episode |
| `po.test` | `po_test` | Translated for 2 to 6 columns |
| `portfolio.optim` | `portfolio_optimize` | Translated; optimizer differs from `quadprog` |
| `pp.test` | `pp_test` | Translated |
| `quadmap` | `quadratic_map` | Translated |
| `runs.test` | `runs_test` | Translated for integer binary vectors |
| `sharpe` | `sharpe_ratio` | Translated using first differences, matching R |
| `sterling` | `sterling_ratio` | Translated |
| `surrogate` | `permutation_surrogate`, `fft_surrogate`, `amplitude_surrogate` | Core generators translated |
| `terasvirta.test` | `terasvirta_test` | Time-series interface translated |
| `tsbootstrap` | `stationary_bootstrap`, `block_bootstrap` | Sample generation translated; callback statistics omitted |
| `white.test` | `white_test` | Time-series interface translated; reproducible seed added |

## Supporting methods represented by result types

| R methods | Fortran replacement |
|---|---|
| `coef.arma`, `residuals.arma`, `fitted.arma`, `vcov.arma`, `summary.arma` | Components of `arma_result` |
| `coef.garch`, `residuals.garch`, `fitted.garch`, `vcov.garch`, `logLik.garch`, `predict.garch` | Components of `garch_result` and `garch_variance` |
| print methods | Ordinary component access and example programs |

## Deliberately omitted

| R function or feature | Reason |
|---|---|
| `plot.arma`, `plot.garch`, `plot.irts`, `plotOHLC`, `seqplot.ts`, line/point methods | Plotting omitted by request |
| `get.hist.quote` | Internet/data-download functionality omitted |
| `irts`, `as.irts`, `approx.irts`, `value`, `[.irts`, `time.irts` | R-specific irregular-time-series class and interpolation plumbing |
| `read.irts`, `write.irts`, `read.ts`, `read.matrix` | R-specific file and class I/O; a small one-column CSV example is provided instead |
| `weekday`, `daysecond`, `is.businessday`, `is.weekend` | Date/time object utilities rather than numerical core |
| `na.remove` methods | Fortran callers are expected to clean missing data explicitly |
| S3 dispatch, formulas, calls, names, attributes, print and summary formatting | R language/runtime behavior |
| Package datasets | The project focuses on computation rather than bundled data |
