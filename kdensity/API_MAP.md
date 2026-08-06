# API map

## Main estimator

| R computation | Fortran API | Notes |
|---|---|---|
| `kdensity(...)` | `fit_kdensity(x, options)` | Returns `type(kdensity_fit)` |
| custom kernel/start | `fit_kdensity_custom(...)` | Uses procedure-pointer specifications |
| callable returned function | `fit%pdf(y)` / `fit%pdf_vector(y)` | Scalar/vector evaluation |
| parametric start evaluation | `fit%start_pdf(y)` | Equivalent to plotting/evaluating the start |
| `coef.kdensity` | `fit%parameters` | Fitted start parameters |
| `logLik.kdensity` | `fit%parametric_loglik` | Correct log likelihood |

## Kernels

`get_kernel` maps all built-in kernel names. The formulas in `R/kernels.R` are
implemented in `src/kdensity_kernels.f90`.

## Starts

`get_start` maps common `univariateML` start names and package aliases.
`supported_starts` returns the compiled catalog. Dynamic R environment updates
are replaced by constructing `type(kd_start)` values with procedure pointers.

## Bandwidths

| R selector | Fortran API |
|---|---|
| `nrd0` | `bandwidth_nrd0` |
| `nrd` | `bandwidth_nrd` |
| `RHE` | `bandwidth_rhe` |
| `JH` | `bandwidth_jh` |
| `HS` | `bandwidth_hs` |
| `ucv` | `bandwidth_ucv` |
| automatic choice | `standard_bandwidth_name` |
| named dispatcher | `select_bandwidth` |

## Helpers

| R helper | Fortran equivalent |
|---|---|
| `get_kernel_start_support` | logic in `fit_kdensity` and `infer_kernel_name` |
| `support_compatible` | `support_compatible` |
| `get_kernel` | `get_kernel` |
| `get_start` | `get_start` |
| `get_bw` | `select_bandwidth` |
| `get_standard_bw` | `standard_bandwidth_name` |
| `kdensity_sq` | normalized square-density calculation inside UCV |
| `compute_hs_bandwidth` | `bandwidth_hs` |

## Omitted R interfaces

`plot`, `lines`, `points`, `$`, `[[`, replacement methods, `update`, `print`,
and `summary` are R object/graphics infrastructure and are not part of the
numerical Fortran library.
