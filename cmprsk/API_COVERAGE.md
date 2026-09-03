# API coverage

Upstream package: **cmprsk 2.2-12** (2024-05-14).

## Computational coverage

| Upstream API/kernel | Fortran API | Coverage |
| --- | --- | --- |
| `cuminc` | `fit_cuminc` | Numerical behavior implemented for complete numeric/integer inputs, including multiple causes/groups and stratified Gray tests. |
| `cinc` | `cumulative_incidence` | Direct structured translation; step coordinates, estimates, and variances match the original native routine in deterministic parity tests. |
| `crstm` / `crst` | `gray_test` | Direct structured translation including stratification, K groups, covariance, chi-square statistic, p-value, and arbitrary `rho`. |
| `timepoints` / `tpoi` | `cuminc_timepoints`, `curve_timepoints`, `timepoint_indices` | Numerical step-function lookup implemented. |
| `crr` | `fit_crr` | Fine-Gray objective, score, Hessian, censoring weighting, Newton/Armijo fitting, robust covariance, residuals, and baseline jumps implemented. |
| `crrfsv` | internal `crr_objective_score_info` | Direct translation; deterministic output matches upstream native routine. |
| `crrf` | internal `crr_objective` | Direct translation. |
| `crrvv` | internal `crr_variance_kernel` | Direct translation of information and sandwich-meat calculations; deterministic output matches upstream. |
| `crrsr` | internal `crr_score_residuals_kernel` | Direct translation; deterministic output matches upstream. |
| `crrfit` | internal `crr_baseline_jumps_kernel` | Direct translation; deterministic output matches upstream. |
| `predict.crr` | `predict_crr` | Cumulative subdistribution-hazard integration and CIF transformation implemented. |
| `summary.crr` numerical calculations | `summarize_crr` | SEs, z tests, normal p-values, exponentiated effects, confidence intervals, and pseudo-LR statistic implemented. |

## Intentional interface omissions

The following are R-specific presentation or orchestration rather than standalone numerical algorithms and are not translated:

- S3 print, subset, summary-object formatting, and plot methods.
- R formula/model-frame/data-frame/factor handling.
- `subset`/`na.action` dispatch and automatic case omission.
- R callback evaluation for `tf`; callers supply the already evaluated time-function matrix.
- Curve names/dimnames and other R list metadata.

These omissions do not remove the package's competing-risks estimators, tests, regression likelihood, covariance calculations, residuals, baseline fit, or predictions.
