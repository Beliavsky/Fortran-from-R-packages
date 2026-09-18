# API coverage

Coverage is based on distinct exported computational R functions in the upstream `NAMESPACE`. Plotting, drawing, printing, formatting, and presentation-only methods are excluded.

Package status: **substantial**  
Mapped computational functions: **8 of 8 (100%)**

| R function | R source | Fortran procedure | Status | Main compatibility note |
|---|---|---|---|---|
| `kde1d` | `upstream/R/kde1d.R` | `kde1d_api::kde1d_fit` | substantial | Numerical estimator translated; R S3/list and ordered-factor construction omitted; direct convolution replaces FFT. |
| `dkde1d` | `upstream/R/kde1d-methods.R` | `kde1d_api::dkde1d` | substantial | Numerical density/mass operation translated; factor metadata omitted. |
| `pkde1d` | `upstream/R/kde1d-methods.R` | `kde1d_api::pkde1d` | substantial | Numerical CDF operation translated; factor metadata omitted. |
| `qkde1d` | `upstream/R/kde1d-methods.R` | `kde1d_api::qkde1d` | substantial | Numerical quantiles translated; ordered-factor reconstruction omitted. |
| `rkde1d` | `upstream/R/kde1d-methods.R` | `kde1d_api::rkde1d` | substantial | Quantile-transform simulation translated; RNG streams/scrambling differ. |
| `equi_jitter` | `upstream/R/jitter.R` | `kde1d_api::equi_jitter` | substantial | Jitter values translated; R factor class/labels omitted. |
| `logLik.kde1d` | `upstream/R/kde1d-methods.R` | `kde1d_api::kde1d_loglik` | substantial | Scalar value translated; `df` is stored as `model%edf` rather than an R attribute. |
| `summary.kde1d` | `upstream/R/kde1d-methods.R` | `kde1d_api::kde1d_summary` | substantial | Same five numerical summary entries; names/printing/invisible return omitted. |

## Computational kernels translated

The Fortran implementation includes:

- conditionally equidistant jittering for integer category codes;
- linear binning with observation weights;
- Gaussian density derivatives;
- modified Sheather-Jones plug-in bandwidth selection used by kde1d;
- local log-constant, log-linear, and log-quadratic fitting;
- cubic interpolation, numerical CDF integration, and quantile inversion;
- one-sided fourth-root and two-sided regularized-probit support transformations;
- Jacobian correction back to the original scale;
- endpoint-tail classification and nonnegative local-linear boundary repair;
- fusion of boundary and transformed bulk estimates;
- discrete mass normalization and discrete CDF/quantile calculations;
- zero-inflated hurdle density, CDF, quantile, likelihood, and all-zero special case;
- weighted log-likelihood and effective degrees of freedom;
- pseudo- and quasi-random quantile-transform simulation.

The coverage fraction records mapped exported computational functions, not exact compatibility with R objects, external RNG implementations, or floating-point evaluation order.
