# API coverage

## Coverage basis

Coverage counts distinct exported computational R functions and explicitly registered S3 methods that perform computation. Plotting, printing, formatting, summaries used for presentation, data-only interfaces, trivial metadata accessors, and generics that merely dispatch are excluded.

The resulting denominator is 19 functions. Fifteen have meaningful Fortran mappings, so the package status is **substantial** and the mapped-function fraction is **15/19 = 0.7894736842105263 (79%)**. This is not a claim of full R compatibility.

## Mapped functions

| R function | Status | Fortran implementation | Main compatibility note |
| --- | --- | --- | --- |
| `contr.earth.response` | substantial | `earth_api:contr_earth_response` | Numeric identity matrix; R dimnames omitted. |
| `earth.default` | substantial | `earth_fit_mod:earth_fit` | Numeric matrices supported; R data-frame/factor, GLM, and CV layers omitted. |
| `earth.formula` | partial | `earth_fit_mod:earth_fit` | Caller must construct the numeric design matrix. |
| `earth.fit` | substantial | `earth_fit_mod:earth_fit` | Core MARS fit and backward GCV pruning; advanced options and alternate pruning methods omitted. |
| `evimp` | substantial | `earth_api:earth_evimp` | Core importance scores; R sorting/trimming/names omitted. |
| `expand.bpairs.default` | substantial | `earth_api:expand_bpairs` | Core count expansion and index semantics; data-frame/sort metadata omitted. |
| `expand.bpairs.formula` | partial | `earth_api:expand_bpairs` | Formula processing omitted. |
| `coef.earth` | substantial | `earth_api:earth_coefficients` | Earth coefficients for all responses; GLM/decomposition presentation omitted. |
| `deviance.earth` | complete | `earth_api:earth_deviance` | Returns earth RSS. |
| `extractAIC.earth` | substantial | `earth_api:earth_extract_aic` | Effective parameter count plus GCV, as upstream intentionally does. |
| `hatvalues.earth` | complete | `earth_api:earth_hatvalues` | Returns selected-regression leverages. |
| `model.matrix.earth` | substantial | `earth_basis:earth_model_matrix` | Numeric predictors only; R formula/factor/naming layer omitted. |
| `predict.earth` | substantial | `earth_api:earth_predict` | Response prediction only; GLM, intervals, `type='terms'`, and varmod omitted. |
| `resid.earth` | partial | `earth_api:earth_residuals` | Ordinary earth residuals only. |
| `residuals.earth` | partial | `earth_api:earth_residuals` | Ordinary earth residuals only. |

## Untranslated computational functions

- `coef.varmod` — transforms coefficients of earth's residual variance submodel; the varmod subsystem is not translated.
- `mars.to.earth` — converts an `mda::mars` R object to an earth S3 object; this is tightly coupled to R object layouts and mda compatibility.
- `predict.varmod` — variance-model standard errors/interval calculations; the varmod subsystem is not translated.
- `update.earth` — reconstructs and reevaluates an R call, formula, data, subset, and pruning arguments; there is no direct Fortran analogue yet.

## Explicitly excluded from the denominator

`earth`, `expand.bpairs` are dispatch-only generics. `anova.earth` and `effects.earth` only warn and return `NULL`. `fitted.earth`, `fitted.values.earth`, `family.earth`, `case.names.earth`, `variable.names.earth`, and `weights.earth` are accessors/delegators rather than distinct numerical kernels. `summary.earth` and `summary.varmod` are presentation-oriented object summaries. All plot/print/format methods and `earth_plotmodsel`, `plot.earth.models`, and `plotd` are presentation-only and intentionally skipped.
