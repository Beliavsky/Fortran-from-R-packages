# dccmidas

Modern free-form Fortran translation of the computational core of the R package **dccmidas 0.1.3** by Vincenzo Candila.

The package provides multivariate volatility/correlation kernels for corrected DCC, asymmetric DCC, DCC-MIDAS, asymmetric DCC-MIDAS, DECO, scalar/diagonal BEKK, RiskMetrics, moving covariance, covariance-loss evaluation, QML standard errors, and substantial in-sample fitting workflows.

## Build

This package is intended to live as a top-level sibling in `Fortran-from-R-packages` with these translated/shared dependencies available:

- `../rfortran-core`
- `../rfortran-linalg`
- `../roll`
- `../rumidas`
- `../rugarch`
- `../maxLik`

Then run:

```text
fpm build
fpm test
fpm run --example basic_dccmidas
```

No separately installed BLAS/LAPACK libraries are requested by this manifest. Linear algebra is reached through `rfortran-linalg`; dependency source is not copied into this package.

## Data conventions

The Fortran interface uses dense arrays instead of R list/xts objects:

- return matrices: `(time, assets)`
- standardized residuals: `(assets, time)`
- conditional covariance/correlation arrays: `(assets, assets, time)`
- `Dt`: `(assets, assets, time)` diagonal standard-deviation matrices

Callers should perform any date alignment or missing-data merging before calling the Fortran routines.

## High-level fitting

`dcc_fit` implements the in-sample two-step workflow. Standard first-stage models reuse the sibling `rugarch` translation for `sGARCH`, `gjrGARCH`, `eGARCH`, `iGARCH`, and `csGARCH`. `norm` and `std` are currently mapped. The `GM_noskew`, `GM_skew`, `DAGM_noskew`, and `DAGM_skew` branches reuse sibling `rumidas`; these branches accept already prepared macro lag matrices with shape `(K+1,time,assets)`.

The second stage uses sibling `maxLik` with the upstream model constraints and observation-wise finite-difference scores for the QML sandwich covariance. Correlation choices are `cDCC`, `aDCC`, `DECO`, `DCCMIDAS`, and `ADCCMIDAS`.

`bekk_fit` fits scalar or diagonal BEKK with sibling `maxLik`. Unlike upstream R, which searches random Uniform starts, the Fortran wrapper uses a deterministic default start unless a start vector is supplied.

Out-of-sample rolling forecast orchestration and R S3/date metadata are intentionally omitted.

## Numerical compatibility notes

The direct model recursions preserve package-specific conventions, including the upstream RiskMetrics lambda placement, raw-outer-product moving covariance, corrected-DCC lagged diagonal scaling, asymmetric intercept adjustments/indexing, and the exact covariance-loss formulas used by `cov_eval`.

The RcppArmadillo `Det` and `Inv` helpers map to `rfortran-linalg`; Rcpp/Armadillo is not vendored or required by this Fortran package.

See `PROVENANCE.md` for detailed compatibility differences and retained upstream sources.

## Translation coverage

Package status: **substantial**.

**22 of 22 (100%)** exported computational R functions are mapped. This fraction measures mapped computational functions, not complete R interface or behavioral compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `dccmidas_loglik` | `dccmidas_dcc` | `dccmidas_loglik` | substantial |
| `dccmidas_mat_est` | `dccmidas_dcc` | `dccmidas_mat_est` | substantial |
| `riskmetrics_mat` | `dccmidas_evaluation` | `riskmetrics_mat` | complete |
| `moving_cov` | `dccmidas_evaluation` | `moving_cov` | complete |
| `sBEKK_loglik` | `dccmidas_bekk` | `sBEKK_loglik` | complete |
| `sBEKK_mat_est` | `dccmidas_bekk` | `sBEKK_mat_est` | complete |
| `dBEKK_loglik` | `dccmidas_bekk` | `dBEKK_loglik` | complete |
| `dBEKK_mat_est` | `dccmidas_bekk` | `dBEKK_mat_est` | complete |
| `dcc_loglik` | `dccmidas_dcc` | `dcc_loglik` | complete |
| `dcc_mat_est` | `dccmidas_dcc` | `dcc_mat_est` | complete |
| `a_dcc_loglik` | `dccmidas_dcc` | `a_dcc_loglik` | complete |
| `a_dcc_mat_est` | `dccmidas_dcc` | `a_dcc_mat_est` | complete |
| `a_dccmidas_loglik` | `dccmidas_dcc` | `a_dccmidas_loglik` | substantial |
| `a_dccmidas_mat_est` | `dccmidas_dcc` | `a_dccmidas_mat_est` | substantial |
| `deco_loglik` | `dccmidas_dcc` | `deco_loglik` | complete |
| `deco_mat_est` | `dccmidas_dcc` | `deco_mat_est` | complete |
| `QMLE_sd` | `dccmidas_evaluation` | `qmle_sd` | complete |
| `dcc_fit` | `dccmidas_fit` | `dcc_fit, dcc_fit_second_stage` | substantial |
| `cov_eval` | `dccmidas_evaluation` | `cov_eval` | complete |
| `bekk_fit` | `dccmidas_fit` | `bekk_fit` | substantial |
| `Det` | `dccmidas_evaluation` | `det_matrix` | complete |
| `Inv` | `dccmidas_evaluation` | `inv_matrix` | complete |

Presentation-only `print.dccmidas`, `summary.dccmidas`, and `plot_dccmidas` are excluded from the computational coverage basis.

## License and provenance

The upstream package declares GPL-3. This translation is distributed under GPL-3.0-only; see `LICENSE`, `NOTICE.md`, and `PROVENANCE.md`. Original authorship and citation material are retained under `upstream/`.
