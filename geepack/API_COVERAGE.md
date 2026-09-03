# API coverage

This document maps the reusable computational surface of R `geepack` 1.3.13 to
the Fortran API. R formula parsing, data-frame/model-frame construction, S3
method dispatch, printing, plotting, and broom/tidy integration are intentionally
out of scope.

| Upstream surface | Fortran coverage | Notes |
| --- | --- | --- |
| `geese.fit` / numerical `geese` engine | **Implemented** | `fit_geese`; simultaneous mean, scale, and association equations |
| `geeglm` numerical fit | **Implemented** | Same GEE engine after caller supplies the numeric design matrix |
| Mean links | **Implemented** | identity, logit, probit, cloglog, log, reciprocal, Fisher-z, LWYBC2, LWYLOG |
| Variance functions | **Implemented** | Gaussian, binomial, Poisson, Gamma |
| Working correlations | **Implemented** | independence, exchangeable, AR(1), unstructured, user-defined, fixed |
| Wave-specific links/variance | **Implemented** | Per-wave integer code arrays in `gee_spec` |
| Mean offsets and weights | **Implemented** | Numeric vectors |
| Scale submodel | **Implemented** | `zsca`, scale links, fixed or estimated scale |
| Association submodel | **Implemented** | `zcor`, association link, fixed/user-defined structures |
| Robust/naive covariance | **Implemented** | beta/alpha/gamma blocks |
| Stable alpha covariance | **Implemented** | Association-only sandwich analogue |
| Cluster influence | **Implemented** | `gee_result%influence` |
| AJS/J1S/FIJ jackknife covariance | **Implemented** | Ordinary GEE through `gee_spec` flags |
| `genZcor` | **Implemented** | `gen_zcor` |
| ordinal `genZodds` | **Implemented** | `gen_zodds` |
| `fixed2Zcor` | **Implemented** | `fixed_to_zcor` |
| `ordgee` | **Implemented** | Cumulative indicators, logit/probit/cloglog, local-odds association, robust/naive covariance |
| QIC / CIC / QICu / QICC arithmetic | **Implemented** | `compute_qic`; caller supplies the independence-model naive beta covariance rather than invoking an R refit wrapper |
| `compCoef` numerical comparison | **Implemented** | `compare_coefficients` uses fitted influences |
| `relRisk` COPY construction/fitting | **Implemented** | `make_relative_risk_copy`, `fit_relative_risk`; final GEE target matches upstream, initialization is explicit/safe rather than R `glm.fit` orchestration |
| coefficient summaries | **Implemented** | Standard errors, Wald chi-square, shared chi-square upper-tail probability |
| `anova.geeglm` numerical Wald kernel | **Implemented** | `wald_contrast`; R nested-formula/model-object orchestration omitted |
| fitted values | **Implemented** | `fitted_means` |
| print/summary formatting | **Excluded** | R presentation layer |
| plotting | **Excluded** | R plotting layer |
| formula/model-frame/NA-action handling | **Excluded** | Caller supplies already prepared numeric arrays |
| broom/tidy methods | **Excluded** | R-specific interface layer |

## Data conventions

- Observations are rows and must be arranged as contiguous clusters.
- `cluster_sizes` gives the number of rows in each cluster and must sum to the
  response length.
- Wave/visit numbers are one-based positive integers used for wave-specific links.
- Optional `cor_param` is real-valued correlation metadata; AR(1) uses its pairwise distances, matching upstream `corp`.
- For unstructured working correlation, the parameter vector uses the global upper-triangle ordering.
- Categorical `fit_ordgee` responses are integer levels `1:nlevels`.
- The public numerical API returns status/error codes instead of throwing R
  conditions.

## Remaining differences

There is no missing major standalone statistical kernel identified in the
upstream package. Remaining differences are interface/orchestration details:
R performs formula parsing, automatic coefficient initialization through GLM
fits, data reordering/NA handling, automatic independence refits inside QIC,
and S3 result formatting. These are intentionally kept outside the numerical
Fortran package.
