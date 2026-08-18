# API map

The upstream package uses `exportPattern("^[[:alpha:]]+")`.

## Computational exports translated

| R function | Fortran mapping |
|---|---|
| `cofi` | `cofi` |
| `dhermite` | `dhermite`; exact path also available as `dhermite_exact` |
| `edg` | `edg` |
| `int.hermite` | `int_hermite` |
| `phermite` | `phermite` |
| `qhermite` | `qhermite` |
| `rhermite` | `rhermite` |
| `glm.hermite` | `fit_glm_hermite` using a numeric design matrix |

The raw regression probability and log-likelihood kernels are additionally
available as `hermite_prob_mu_d` and `hermite_glm_loglik`.

## R-interface/presentation exports omitted

- `print.summary.glm.hermite`
- `summary.glm.hermite`

The Fortran fit result already contains the numerical values needed to build
such a summary: coefficients, covariance/Hessian, fitted means, log-likelihood,
likelihood-ratio statistic, p-value, and AIC. Formula parsing, model frames,
factor contrasts, and S3 printing are intentionally not emulated.
