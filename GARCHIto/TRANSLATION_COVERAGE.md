# Translation coverage

Upstream `NAMESPACE` exports three routines. All three are represented.

| Upstream export | Status | Fortran API |
|---|---|---|
| `UnifiedEst` | translated | `unified_est` |
| `RealizedEst` | translated | `realized_est` |
| `RealizedEst_Option` | translated | `realized_est_option` |

## Internal computational code

The following nested R computations are included in typed Fortran form:

- unified conditional variance recursion and quasi-likelihood;
- realized no-jump recursion and quasi-likelihood;
- realized jump recursion, median-jump initialization, and quasi-likelihood;
- homogeneous and heterogeneous option measurement likelihoods;
- stationarity and box constraints;
- OLS initialization of `a`, `b`, and `sigma_e`;
- fitted conditional variance and one-step forecast construction.

## Intentionally not translated

- R package loading and `.rda` object infrastructure;
- roxygen/S3-style documentation mechanics;
- vignette plotting and HTML generation;
- the external Rsolnp implementation itself.

No exported plotting routine or other exported computational routine is
omitted.
