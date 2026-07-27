# Computational coverage

This file maps the exported numerical surface of SharpeR 1.4.0 to the modern
Fortran implementation. R class adapters and display methods are described
separately because they do not represent distinct numerical algorithms.

## Sharpe-ratio objects and estimation

| SharpeR interface | Fortran interface | Status |
| --- | --- | --- |
| `sr`, `as.sr` for numeric vectors/matrices | `fit_sr`, `fit_sr_matrix` | Implemented |
| `is.sr` | Fortran type checking through `type(sr_result)` | Native equivalent |
| `reannualize.sr` | `reannualize_sr` | Implemented |
| `se.sr` | `sr_standard_error` | Implemented |
| `confint.sr` | `sr_confint` | Implemented |
| `predint.sr` | `predint` | Implemented |
| `sr_bias` | `sr_bias` | Implemented |
| `sr_variance` | `sr_variance` | Implemented |
| `sr_vcov` | `sr_vcov` | Implemented |

The higher-order paths include k-statistics and standardized cumulants used by
SharpeR's nonnormal bias and variance approximations.

## Optimal Sharpe ratios and spanning

| SharpeR interface | Fortran interface | Status |
| --- | --- | --- |
| `sropt`, `as.sropt` | `fit_sropt`, `make_sropt` | Implemented |
| `is.sropt` | `type(sropt_result)` | Native equivalent |
| `reannualize.sropt` | `reannualize_sropt` | Implemented |
| `confint.sropt` | `sropt_confint` | Implemented |
| `asnr_confint.sropt` | `achieved_snr_confint` | Implemented |
| `inference.sropt` | `infer_sropt` | Implemented |
| `sric` | `sric` | Implemented |
| `del_sropt`, `as.del_sropt` | `fit_del_sropt`, `make_del_sropt` | Implemented |
| `is.del_sropt` | `type(del_sropt_result)` | Native equivalent |
| `confint.del_sropt` | `sropt_confint` on the delta representation | Implemented |
| `asnr_confint.del_sropt` | `achieved_delta_snr_confint` | Implemented |
| `inference.del_sropt` | `infer_del_sropt` | Implemented |

Markowitz weights and the full/subspace/delta Hotelling T-squared quantities
are returned directly in typed result objects.

## Distribution functions

| SharpeR interface | Fortran interface | Status |
| --- | --- | --- |
| `dsr`, `psr`, `qsr`, `rsr` | Same names | Implemented |
| `dsropt`, `psropt`, `qsropt`, `rsropt` | Same names | Implemented |
| `plambdap`, `qlambdap`, `rlambdap` | Same names | Implemented |
| `pco_sropt`, `qco_sropt` | Same names | Implemented |

The internal rescaled-t, T-squared, noncentral t, noncentral F, noncentral
chi-square, normal, Student-t, chi-square, and F routines are also included.

## Tests and power

| SharpeR interface | Fortran interface | Status |
| --- | --- | --- |
| `sr_test` | `sr_test` | Implemented |
| paired Sharpe comparison path | `paired_sr_test` | Implemented |
| `sr_unpaired_test` | `unpaired_sr_test` | Implemented |
| `sr_equality_test` | `sr_equality_test` | Implemented |
| `sropt_test` | `sropt_test` | Implemented |
| `sr_max_test` | `sr_max_test` | Implemented |
| `sr_conditional_test` | `sr_conditional_test` | Implemented with native polyhedral inputs |
| `power.sr_test` | `power_sr_test`, `required_n_sr_test` | Implemented |
| `power.sropt_test` | `power_sropt_test`, `required_df2_sropt_test` | Implemented |

## Unified moment covariance

| SharpeR interface | Fortran interface | Status |
| --- | --- | --- |
| `sm_vcov` | `sm_vcov` | Implemented |
| `ism_vcov` | `ism_vcov` | Implemented |

The empirical route uses sample covariances of the stacked first and second
moments. The normal-model route uses the corresponding Gaussian analytic
covariance formulas. `ism_vcov` applies a numerical delta-method Jacobian to
the inverse-second-moment transformation.

## R-specific infrastructure

The following are not compiled because they are class, presentation, or
external-environment adapters rather than independent numerical procedures:

- S3 `print` and `summary` methods.
- `xts`, `zoo`, `timeSeries`, data-frame, and `lm` adapters.
- R formula and missing-value dispatch.
- Package data objects and vignettes.
- R table or plotting helpers.

Equivalent calculations can be performed by passing ordinary Fortran arrays
and explicit scalar options to the typed APIs.
