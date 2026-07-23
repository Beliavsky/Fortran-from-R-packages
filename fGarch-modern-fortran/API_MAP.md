# API map

| R fGarch API | Fortran API | Status |
|---|---|---|
| `dstd`, `pstd`, `qstd`, `rstd` | same names | implemented |
| `dged`, `pged`, `qged`, `rged` | same names | implemented |
| `dsnorm`, `psnorm`, `qsnorm`, `rsnorm` | same names | implemented |
| `dsstd`, `psstd`, `qsstd`, `rsstd` | same names | implemented |
| `dsged`, `psged`, `qsged`, `rsged` | same names | implemented |
| base normal helpers | `dnorm_fg`, `pnorm_fg`, `qnorm_fg`, `rnorm_fg` | implemented |
| `absMoments` | `absolute_moment` | analytic supported families |
| `garchSpec` | `garch_spec`, `make_garch_spec` | typed replacement |
| `garchSim` | `simulate_garch` | implemented; arbitrary p/q and optional ARMA terms |
| `garchFit` | `fit_garch11`, `fit_aparch11` | partial: 1,1 fitters only |
| internal likelihood | `garch_log_likelihood`, `garch_filter` | implemented |
| `garchKappa` | `garch_kappa` | numerical integration |
| `.truePersistence` | `true_persistence` | implemented |
| `predict` volatility | `forecast_volatility` | APARCH/GARCH forecast |
| `VaR` | `value_at_risk` | implemented |
| `ES` | `expected_shortfall` | implemented numerically |
| distribution fit functions | `fit_distribution` | common MLE replacement |
| `tsdiag` pieces | `jarque_bera_statistic`, `ljung_box_statistic` | partial |
| `garchFitControl` | procedure arguments | not represented as an object |
| optimizer and Hessian methods | internal Nelder-Mead | partial; Hessian omitted |
| `snig` conditional distribution | none | omitted |
| `.gogarchFit` / GO-GARCH | none | omitted |
| plotting methods and sliders | none | omitted as requested |
| S3/S4 methods, formulas, update/show/summary | derived types and examples | R-specific code omitted |
| timeSeries and bundled data integration | none | omitted |

## Public modules

Most users should write:

```fortran
use fgarch
```

The umbrella module re-exports the public procedures and types from the component modules.
