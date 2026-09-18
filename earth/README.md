# earth — modern Fortran translation

This directory translates the computational core of the R package **earth 5.3.6** to modern free-form Fortran with FPM packaging. The main target is multivariate adaptive regression splines (MARS): forward basis construction, weighted least-squares fitting, backward pruning with earth's GCV definition, prediction/model matrices, regression diagnostics, binomial-pair expansion, and variable importance.

The public Fortran API is in `earth_api`. All maintained real arithmetic uses the single `dp = real64` kind from `earth_kinds`. The implementation has no BLAS/LAPACK dependency; its least-squares path is self-contained, so no system numerical libraries or vendored shared dependencies are required.

## Build and test

With FPM installed:

```text
fpm build
fpm test
fpm run --example basic_earth
```

The intended deployment is a top-level `earth/` directory in `Beliavsky/Fortran-from-R-packages`.

## Main Fortran API

- `earth_fit` — numeric-matrix MARS fitting with degree, penalty, basis-count/threshold controls, automatic or explicit minspan/endspan, new-variable penalty, forced linear predictors, case weights, multiple responses, backward GCV pruning, and no-prune selection.
- `earth_predict` and `earth_model_matrix` — evaluate the selected regression model on new numeric predictor matrices.
- `earth_coefficients`, `earth_residuals`, `earth_hatvalues`, `earth_deviance`, `earth_extract_aic` — selected-model diagnostics and accessors.
- `earth_evimp` — `nsubsets`, GCV, and RSS variable importance from the backward pruning path.
- `expand_bpairs` — expand two-column binomial-pair counts to long form, including earth's zero-zero rule.
- `contr_earth_response` — identity response contrasts.

`earth_model` stores forward directions and cuts, selected term indices, pruning path, coefficients, fitted values, residuals, leverages, RSS/GCV paths, and fit summary fields.

## Numerical scope and differences from R earth

The forward pass follows earth/Friedman MARS conventions for hinge pairs, hierarchical interactions, linear-predictor restrictions, span heuristics, GCV penalties, and new-variable penalties, but uses direct candidate regressions rather than the upstream C implementation's incremental covariance/cache machinery. Rank-deficient hinge pairs can retain one useful child, corresponding to the practical effect of earth's regression-fix path.

The current translation intentionally does **not** reproduce R formula/model-frame/factor handling, GLM post-processing, cross-validation/bagging, variance models and prediction intervals, response weights `wp`, `allowed` callbacks, `fast.k`/`fast.beta` acceleration controls, `Auto.linpreds`, or the non-backward pruning methods. R S3 object construction, warning text, names/attributes, plotting, printing, and other presentation behavior are also outside scope. The numeric least-squares implementation uses reorthogonalized modified Gram-Schmidt instead of upstream `leaps`/AS274 code.

## Translation coverage

Package status: **substantial**. Coverage is **15 of 19 (79%)** exported computational R functions.

The fraction measures R functions with a meaningful Fortran mapping, **not** complete R-interface compatibility. Plotting, printing, formatting, presentation-only methods, simple accessors, and dispatch-only generics are excluded from the denominator under the project's coverage convention.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `contr.earth.response` | `earth_api` | `contr_earth_response` | substantial |
| `earth.default` | `earth_fit_mod` | `earth_fit` | substantial |
| `earth.formula` | `earth_fit_mod` | `earth_fit` | partial |
| `earth.fit` | `earth_fit_mod` | `earth_fit` | substantial |
| `evimp` | `earth_api` | `earth_evimp` | substantial |
| `expand.bpairs.default` | `earth_api` | `expand_bpairs` | substantial |
| `expand.bpairs.formula` | `earth_api` | `expand_bpairs` | partial |
| `coef.earth` | `earth_api` | `earth_coefficients` | substantial |
| `deviance.earth` | `earth_api` | `earth_deviance` | complete |
| `extractAIC.earth` | `earth_api` | `earth_extract_aic` | substantial |
| `hatvalues.earth` | `earth_api` | `earth_hatvalues` | complete |
| `model.matrix.earth` | `earth_basis` | `earth_model_matrix` | substantial |
| `predict.earth` | `earth_api` | `earth_predict` | substantial |
| `resid.earth` | `earth_api` | `earth_residuals` | partial |
| `residuals.earth` | `earth_api` | `earth_residuals` | partial |

Untranslated computational functions in this coverage basis are `coef.varmod`, `mars.to.earth`, `predict.varmod`, and `update.earth`. See `API_COVERAGE.md` for the counting rationale and `fpm.toml` for machine-readable mappings.

## Provenance and licensing

The upstream R package metadata and source snapshot used for this translation are retained under `upstream/` for provenance and are not compiled by FPM. See `NOTICE.md`, `PROVENANCE.md`, and `COPYING` for authorship, source history, and license details.
