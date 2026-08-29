# API mapping

| R / upstream computational entry point | Fortran entry point | Notes |
|---|---|---|
| `dtmvnorm` | `dtmvnorm`, `dtmvnorm_one` | Matrix/scalar forms; optional log density |
| `ptmvnorm` | `ptmvnorm` | Probability inside an inner rectangle conditional on truncation |
| `dtmvnorm.marginal` | `dtmvnorm_marginal` | One-dimensional marginal density |
| `dtmvnorm.marginal2` | `dtmvnorm_marginal2` | Bivariate marginal density |
| `ptmvnorm.marginal` | `ptmvnorm_marginal` | One-dimensional marginal CDF |
| `qtmvnorm.marginal` | `qtmvnorm_marginal` | Lower/upper/both-tail bisection |
| `mtmvnorm` | `mtmvnorm` | `tmvnorm_moments_t` result |
| `JohnsonKotzFormula` | internal branch of `mtmvnorm` | Automatically used for partial truncation |
| `mtmvnorm.quadrature` | not separately exposed | Main analytical moment algorithm is translated; diagnostic R `integrate()` wrapper omitted |
| `rtmvnorm` | `rtmvnorm` | High-level rejection/Gibbs wrapper |
| `rtmvnorm.rejection` | `rtmvnorm_rejection` | Supports rectangular or supplied `D` constraints |
| `rtnorm.gibbs` | `rtnorm` | Direct inverse-CDF truncated normal generation |
| `rtmvnorm.gibbs` | `rtmvnorm_gibbs` | Dense covariance form |
| `rtmvnorm.gibbs.Precision` | `rtmvnorm_gibbs_precision` | Dense precision form |
| native `rtmvnormgibbscov` | `rtmvnorm_gibbs` | Modern typed implementation |
| native `rtmvnormgibbsprec` | `rtmvnorm_gibbs_precision` | Modern typed implementation |
| `rtmvnorm2` / `rtmvnorm.gibbs2` | `rtmvnorm2`, `rtmvnorm_gibbs_linear` | General `r x d` linear constraints, including `r>d` |
| native `rtmvnormgibbscov2` / `rtmvnormgibbsprec2` | `rtmvnorm_gibbs_linear` | Common precision-conditional implementation |
| `rtmvnorm.sparseMatrix` / native sparse routines | `rtmvnorm_sparse_csc`, `rtmvnorm_sparse_triplet` | Native array APIs replace Matrix-class wrappers |
| `dtmvt` | `dtmvt`, `dtmvt_one` | Shifted multivariate t parameterization |
| `ptmvt` | `ptmvt` | Conditional rectangle probability |
| `ptmvt.marginal` | `ptmvt_marginal` | One-dimensional marginal CDF |
| `rtmvt` | `rtmvt` | High-level rejection/Gibbs wrapper |
| `rtmvt.rejection` | `rtmvt_rejection` | Rejection sampler |
| `rtmvt.gibbs` | `rtmvt_gibbs` | Geweke-style Gibbs construction |
| `mle.tmvnorm` | `mle_tmvnorm` | Typed result object; optional Cholesky parameterization |
| `gmultiManjunathWilhelm` | `gmm_moments_manjunath_wilhelm` | Returns observation-by-moment matrix |
| `gmultiLee` | `gmm_moments_lee` | Configurable `lmax` |
| `gmm.tmvnorm` | `gmm_tmvnorm` | Two-step GMM with IID moment covariance |

## Interface substitutions

R named `fixed=` parameter lists are represented by an optional logical `free_mask` over the packed full parameter vector. This keeps the computational ability to hold arbitrary parameters fixed without recreating R's dynamic argument-name system.

The sparse APIs use ordinary integer/value arrays rather than R `Matrix` objects. `rtmvnorm_sparse_csc` expects **1-based** row indices and **1-based** column pointers, consistent with normal Fortran indexing. The original R wrapper's zero-based `dgCMatrix` storage conversion is therefore not reproduced.

`gmm.tmvnorm` upstream delegates arbitrary `...` options to the external R package `gmm`. The Fortran port translates tmvtnorm's two moment systems and provides a self-contained conventional two-step IID GMM estimator. The lower-level moment matrices are public so callers can apply alternative HAC/weighting machinery if desired.
