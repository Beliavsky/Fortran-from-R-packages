# grpreg — modern Fortran translation

This package translates the computational core of Patrick Breheny and collaborators' R package **grpreg 3.6.0** to modern free-form Fortran with FPM.

It provides grouped regularization paths for Gaussian regression, binomial and Poisson GLMs, and Cox proportional-hazards models. Implemented penalties include group lasso (`grLasso`), group MCP (`grMCP`), group SCAD (`grSCAD`), group exponential lasso (`gel`), composite MCP (`cMCP`), and group bridge (`gBridge`). The implementation is self-contained and does not require BLAS/LAPACK or R at build time.

## Build

```text
fpm build
fpm test
fpm run --example basic_usage
```

The public umbrella module is `grpreg`; the package-wide real kind `dp` is re-exported from it.

A minimal fit looks like:

```fortran
use grpreg

type(grpreg_fit_type) :: fit
call grpreg_fit(x, y, group, fit, penalty='grMCP', nlambda=100)
```

The numerical API also includes group bridge, grouped Cox fitting and survival prediction, CV, mFDR, information-criterion selection, residuals/log likelihood, spline expansion, and the nonlinear demonstration-data generator.

## Numerical design

Predictors are centered and scaled with the upstream population-norm convention. Group-selection penalties additionally orthogonalize each group so `X_g^T X_g / n = I`, matching the condition used by the upstream closed-form group updates. Rank-deficient group directions are removed. Coefficients are transformed back to the original predictor scale before being returned.

Strong rules/BEDPP screening are deliberately not required: every group is checked in the descent loop. This changes speed, not the stated penalized objective. The package contains a small symmetric Jacobi eigensolver rather than depending on BLAS/LAPACK.

## Translation coverage

**Package status: substantial. Coverage: 19 of 19 (100%).**

The percentage measures mapped exported computational functions/methods, **not complete R compatibility**.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `grpreg` | `grpreg_api` | `grpreg_fit` | substantial |
| `grpsurv` | `grpreg_api` | `grpsurv_fit` | substantial |
| `gBridge` | `grpreg_api` | `gbridge_fit` | substantial |
| `cv.grpreg` | `grpreg_cv` | `cv_grpreg` | substantial |
| `cv.grpsurv` | `grpreg_cv` | `cv_grpsurv` | substantial |
| `expand_spline` | `grpreg_spline` | `expand_spline`, `predict_spline` | substantial |
| `gen_nonlinear_data` | `grpreg_data` | `gen_nonlinear_data` | substantial |
| `mfdr` | `grpreg_mfdr` | `mfdr_grpreg` | substantial |
| `AUC.cv.grpsurv` | `grpreg_cv` | `auc_cv_grpsurv` | partial |
| `coef.grpreg` | `grpreg_api` | `coef_grpreg` | complete |
| `coef.cv.grpreg` | `grpreg_cv` | `coef_cv_grpreg` | complete |
| `logLik.grpreg` | `grpreg_api` | `loglik_grpreg` | substantial |
| `logLik.grpsurv` | `grpreg_api` | `loglik_grpreg` | complete |
| `predict.grpreg` | `grpreg_api` | `predict_grpreg`, `coef_grpreg`, `count_nonzero`, `count_nonzero_groups`, `group_norms` | substantial |
| `predict.cv.grpreg` | `grpreg_cv` | `predict_cv_grpreg` | substantial |
| `predict.grpsurv` | `grpreg_api` | `predict_grpsurv_link`, `predict_grpsurv_survival`, `predict_grpsurv_hazard`, `predict_grpsurv_median` | substantial |
| `residuals.grpreg` | `grpreg_api` | `residuals_grpreg` | substantial |
| `select.grpreg` | `grpreg_api` | `select_grpreg` | substantial |
| `summary.cv.grpreg` | `grpreg_cv` | `summarize_cv_grpreg` | substantial |

See `API_COVERAGE.md` and the machine-readable `[extra.translation]` tables in `fpm.toml`.

## Material compatibility differences

- Multitask/seemingly-unrelated regression via matrix-valued `y` is not yet translated; the Fortran fit API currently accepts a response vector.
- Upstream strong-rule/BEDPP screening is omitted. The unscreened descent targets the same objective but can be slower on very wide problems.
- `cv.grpsurv` implements the quick SE calculation; the optional bootstrap SE path is not included.
- Default CV folds are deterministic. User-provided folds are supported. This is an RNG/interface difference, not a model difference.
- Natural cubic spline columns span the same boundary-constrained spline space as `splines::ns`, but may be an orthonormal rotation of R's particular basis. B-spline construction follows the standard `bs` basis and is independently parity-tested.
- `AUC.cv.grpsurv` uses explicit comparable-pair concordance. Fine details of `survival::concordancefit` tie weighting are not reproduced.
- `select.grpreg(..., smooth=TRUE)` is not reproduced; all five unsmoothed information criteria are implemented.
- R formula/model-frame handling, S3 classes, names/dimnames, printing, plotting, and callback-style interfaces are intentionally outside the Fortran API.
- Fortran RNG streams do not match R's RNG for `gen_nonlinear_data`.

## Validation

`test/test_grpreg.f90` exercises every translated model family and the package-level helpers. `tools/parity_driver.f90` plus `tools/check_parity.py` provide independent NumPy/SciPy checks for a grouped Gaussian fit, binomial and Poisson likelihood fits, Cox partial likelihood, and B-spline construction.

See `VALIDATION.md` for the exact compiler commands and measured differences from the final release archive.

## License and provenance

The upstream package is GPL-3. The full license text is in `LICENSE`; upstream metadata and source material used during the translation are under `upstream/`. See `NOTICE.md` for attribution and provenance details.
