# statmod R-to-Fortran API mapping

The Fortran API is numerical rather than an emulation of R formula/S3 objects.
Names use lowercase snake_case.

| Upstream R entry point | Fortran computational entry point | Notes |
|---|---|---|
| `logmdigamma` | `logmdigamma` | Stable asymptotic/recurrence implementation. |
| `cumulant.digamma` | `cumulant_digamma` | |
| `meanval.digamma` | `meanval_digamma` | |
| `d2cumulant.digamma` | `d2cumulant_digamma` | |
| `canonic.digamma` | `canonic_digamma` | Three Newton refinements as upstream. |
| `varfun.digamma` | `varfun_digamma` | |
| `unitdeviance.digamma` | `unitdeviance_digamma` | |
| `Digamma` | above numerical kernels | R family-object construction omitted. |
| `dinvgauss` | `dinvgauss` | |
| `pinvgauss` | `pinvgauss` | Stable lower/upper-tail formulas. |
| `qinvgauss` | `qinvgauss` | Newton inversion and upstream starting regimes. |
| `rinvgauss` | `rinvgauss` | Uses `r_mod` RNG. |
| `gauss.quad` | `gauss_quad` | Legendre, Chebyshev I/II, Hermite, Laguerre, Jacobi. |
| `gauss.quad.prob` | `gauss_quad_prob` | Uniform, normal, beta and gamma probability measures. |
| `expectedDeviance` | `expected_deviance`; family-specific routines | Binomial is exposed as `expected_deviance_binomial(prob,nsize,...)`. |
| `glmgam.fit` | `glmgam_fit` | Identity-link Gamma secure fit. |
| `glmnb.fit` | `glmnb_fit` | Log-link negative-binomial secure fit. |
| `fitNBP` | `fit_nbp` | Multi-group common-overdispersion fit. |
| `mixedModel2Fit`, `randomizedBlockFit` | `mixed_model2_fit` | Numerical design-matrix API. |
| `mixedModel2`, `randomizedBlock` | `mixed_model2_fit` | Formula/model-frame construction omitted. |
| `remlscore` | `remlscore` | Mean/dispersion REML scoring. |
| `remlscoregamma` | `remlscoregamma` | Default upstream log-mean/log-dispersion numerical model. |
| `glm.scoretest` | `glm_scoretest` | Accepts residuals, weights, existing and new design matrices. |
| `tweedie` | `tweedie_variance`, `tweedie_linkfun`, `tweedie_linkinv`, `tweedie_mu_eta`, `tweedie_deviance_residual` | R GLM-family object omitted. |
| `eldaOneGroup` | `elda_one_group` | |
| `elda`, `limdil` | `elda_fit` | Integer group labels; formula/print/plot object layer omitted. |
| `matvec` | `matvec` | Column scaling. |
| `vecmat` | `vecmat` | Row scaling. |
| `forward` | `forward_select` | |
| `mscale` | `mscale` | Hampel/Yohai robust scale; helper rho/psi also public. |
| `hommel.test` | `hommel_test` | Literal upstream logical-decision algorithm. |
| `permp` | `permp` | Exact and Gaussian-quadrature approximate paths. |
| `sage.test` | `sage_test` | Exact binomial / large-count chi-square branch. |
| `power.fisher.test` | `power_fisher_test` | Two-sided, less and greater alternatives; uses `r_mod` probability/test helpers and RNG. |
| `meanT` | `mean_t` | |
| `compareTwoGrowthCurves` | `compare_two_growth_curves` | |
| `compareGrowthCurves` | `compare_growth_curves` | Integer groups; Holm adjustment included. |
| `qres.binom` | `qres_binom` | Model-object extraction replaced by explicit numeric arguments. |
| `qres.pois` | `qres_pois` | |
| `qres.gamma` | `qres_gamma` | |
| `qres.invgauss` | `qres_invgauss` | |
| `qres.nbinom` | `qres_nbinom` | Real negative-binomial size supported via beta CDF. |
| `qres.tweedie` | `qres_tweedie` | Uses vendored translated `tweedie` dependency. |
| `qres.default`, `qresid`, `qresiduals` | `qres_default` plus family-specific routines | R S3 dispatch omitted. |
| `plotGrowthCurves`, `plot.limdil`, `print.limdil` | omitted | Presentation/UI code. |
