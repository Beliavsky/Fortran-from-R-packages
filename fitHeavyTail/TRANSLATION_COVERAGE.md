# Translation coverage

## Exported R routines

| Original routine | Fortran routine | Status |
|---|---|---|
| `fit_Tyler` | `fit_tyler` | Translated |
| `fit_Cauchy` | `fit_cauchy` | Translated |
| `fit_mvt` | `fit_mvt` | Translated |
| `fit_mvst` | `fit_mvst` | Translated |
| `nu_OPP_estimator` | `nu_opp_estimator` | Translated |
| `nu_POP_estimator` | `nu_pop_estimator` | Translated |

## Computational helper coverage

The following non-exported numerical functionality is also represented:

- spatial-median initialization and scatter-scale recovery;
- weighted Tyler and Cauchy updates;
- Student-t complete-data and missing-data E-steps;
- factor-analysis loading and idiosyncratic-variance updates;
- multivariate-t and skewed-t log likelihoods;
- ECM/ECME scalar optimization;
- all `nu_mle` modes used in the R source;
- marginal kurtosis, cross cumulants, combined cumulants, Hill, and Pareto
  tail-index estimators;
- skewness estimation;
- generalized-hyperbolic skewed-t conditional expectations;
- real-order Bessel K logarithms, ratios, and order derivatives;
- POP approximate, exact, and sigma-corrected variants;
- deterministic Student-t simulation used by the simulation correction.

## Deliberately omitted

- `plot_convergence` and all `ggplot2`/`reshape2` code;
- roxygen, CRAN, vignette, HTML, and R object infrastructure;
- saved `.RData` regression fixtures;
- iteration-history storage whose sole package use is convergence plotting.

The original R sources and manual pages are retained under `original/` for
traceability.
