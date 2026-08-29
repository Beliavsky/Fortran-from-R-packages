# API mapping

This file maps the main computational pieces of R package `gmm` 1.9-1 to the Fortran numerical API. R formula/model-frame construction, S3 dispatch/presentation, plotting, and printing are intentionally excluded.

| Upstream computational interface | Fortran counterpart |
|---|---|
| `gmm`, `momentEstim.baseGmm` | `gmm_fit`, `gmm_evaluate`, `gmm_fit_fixed_weight` |
| twoStep / iterative / cue | `GMM_TWO_STEP`, `GMM_ITERATIVE`, `GMM_CUE` |
| `tsls` | `tsls_fit` |
| linear GMM internals | `linear_gmm_fit`, `linear_moments`, `linear_gradient` |
| `gel`, `evalGel`, GEL moment estimation | `gel_fit`, `gel_evaluate` |
| `.rho`, `.getCgelLam`, `.Wu`, `.CUE_lam`, RCUE machinery | `gel_rho`, `gel_lambda`, `gel_objective` |
| `.getImpProb` | `gel_implied_prob` |
| EL / ET / CUE / ETEL / HD / ETHD / RCUE | `GEL_EL`, `GEL_ET`, `GEL_CUE`, `GEL_ETEL`, `GEL_HD`, `GEL_ETHD`, `GEL_RCUE` |
| HAC/weight calculations | `moment_covariance`, `hac_covariance`, `kernel_weight`, `smooth_moments` |
| Andrews bandwidth path | `bw_andrews` |
| `bwWilhelm` | `bw_wilhelm` |
| `specTest.gmm`, `specTest.gel` numerical J test | `j_test`; result fields in GMM/GEL fit types |
| `.BigCov`, `KTest` | `kleibergen_k`, `kleibergen_k_from_blocks` |
| `sysGmm` numerical estimator | `system_gmm_fit` |
| `sur` | `sur_fit` |
| `threeSLS` | `three_sls_fit` |
| `five` | `five_fit` |
| random-effect numerical fit | `random_effect_fit` |
| `ategel`, `.momentFctATE`, `.DmomentFctATE` | `ategel_fit`, `ate_moments`, `ate_gradient` |
| ATE marginal effects | `ate_marginal_effects` |
| `charStable` | `char_stable` |
| internal matrix helpers | `gmm_linalg` module |

## Intentionally omitted R interfaces

The following are not numerical algorithms and are not recreated as Fortran object systems: formulas, model frames, `getDat`, S3 `print`/`summary`/`coef`/`residuals` methods, plotting, kernel-object constructors, R time-series attributes, and console/progress behavior.

The Fortran API accepts explicit numeric arrays, integer method selectors, callback procedures, and result derived types instead.
