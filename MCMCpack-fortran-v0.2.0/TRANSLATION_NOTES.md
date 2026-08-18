# Translation notes

## Upstream

- Package: MCMCpack
- Upstream version translated: 1.7-1
- Upstream date: 2024-08-27
- Upstream license: GPL-3
- Primary authors: Andrew D. Martin, Kevin M. Quinn, Jong Hee Park, with additional contributors listed in the original package metadata

Copies of upstream `DESCRIPTION`, `CITATION`, and `README` are retained in `original/`. GPL-3 license text and translation/dependency license material are retained in the release tree.

## Translation scope

The project translates reusable numerical code and sampler engines into modern Fortran. R-only interfaces such as formulas, model frames, factors, S3 classes, printing, plotting, `...` dispatch, `.C`/`.Call` wrappers, R interrupt handling, and R-specific serialization are not reproduced. Fortran entry points therefore accept explicit vectors, matrices, priors, starting values, and control parameters.

### Core and previously translated models

| Upstream area | Fortran entry point/module |
|---|---|
| MCMCregress | `mcmc_regress` |
| MCMCprobit | `mcmc_probit` |
| MCMClogit | `mcmc_logit` |
| MCMCpoisson | `mcmc_poisson` |
| MCMCtobit | `mcmc_tobit` |
| MCMCquantreg | `mcmc_quantreg` |
| MCMCmetrop1R | `mcmc_metrop1r` |
| MCMCmnl | `mcmc_mnl` |
| MCMCnegbin | `mcmc_negbin` |
| MCMCbinaryChange | `mcmc_binary_change` |
| MCMCprobitChange | `mcmc_probit_change` |
| MCMCregressChange | `mcmc_regress_change` |
| MCMCpoissonChange | `mcmc_poisson_change` |
| MCMCnegbinChange | `mcmc_negbin_change` |
| MCMCoprobitChange | `mcmc_oprobit_change` |
| MCMCresidualBreakAnalysis | `mcmc_residual_break_analysis` |
| MCMChregress | `mcmc_hregress` |
| MCMChlogit | `mcmc_hlogit` |
| MCMChpoisson | `mcmc_hpoisson` |
| MCMCirt1d | `mcmc_irt1d` |
| MCMCirtHier1d | `mcmc_irt_hier1d` |
| MCMCirtKd | `mcmc_irtkd` via ordinal-factor kernel |
| MCMCpaircompare | `mcmc_paircompare` |
| MCMCpaircompare2d | `mcmc_paircompare2d` |
| MCMCfactanal | `mcmc_factanal` |
| MCMCordfactanal | `mcmc_ordfactanal` |
| MCMCmixfactanal | `mcmc_mixfactanal` |
| MCMCoprobit | `mcmc_oprobit` |
| MCMCSVDreg | `mcmc_svdreg` |
| MCbinomialbeta | `mc_binomial_beta` |
| MCpoissongamma | `mc_poisson_gamma` |
| MCnormalnormal | `mc_normal_normal` |
| MCmultinomdirichlet | `mc_multinom_dirichlet` |
| dinvgamma/rinvgamma | same names |
| ddirichlet/rdirichlet | same names |
| dwish/rwish, diwish/riwish | same names |
| d/r noncentral hypergeometric | `dnoncenhypergeom`, `rnoncenhypergeom`, `noncenhypergeom_pmf` |
| BayesFactor arithmetic | `bayes_factor` |
| PostProbMod | `post_prob_mod` |
| make.breaklist | `make_breaklist` |
| mptable/topmodels calculations | `marginal_inclusion`, `top_models` |
| vech/xpnd | `vech`, `xpnd` |
| procrustes | `procrustes` |
| WAIC internal utility | `waic` |

### Specialized samplers added in v0.2.0

All sampler engines that were explicitly listed as missing from v0.1.0 now have Fortran implementations:

| Upstream sampler | Fortran entry point | Main module |
|---|---|---|
| MCMChierEI | `mcmc_hier_ei` | `mcmcpack_ei` |
| MCMCdynamicEI | `mcmc_dynamic_ei` | `mcmcpack_ei` |
| MCMCdynamicIRT1d | `mcmc_dynamic_irt1d` | `mcmcpack_dynamic_irt` |
| MCMCirtKdRob | `mcmc_irtkd_rob` | `mcmcpack_irt_robust` |
| MCMCpaircompare2dDP | `mcmc_paircompare2d_dp` | `mcmcpack_paircompare2d_dp` |
| SSVSquantreg | `ssvs_quantreg` | `mcmcpack_ssvs_quantreg` |
| HMMpanelFE | `hmm_panel_fe` | `mcmcpack_panel_hmm` |
| HMMpanelRE | `hmm_panel_re` | `mcmcpack_panel_hmm` |
| HDPHMMpoisson | `hdphmm_poisson` | `mcmcpack_hdp_hmm` |
| HDPHMMnegbin | `hdphmm_negbin` | `mcmcpack_hdp_hmm` |
| HDPHSMMnegbin | `hdphsmm_negbin` | `mcmcpack_hdp_hmm` |
| MCMCirtHier1d PX path | `mcmc_irt_hier1d(..., px=.true.)` | `mcmcpack_irthier` |
| MCMCirtHier1d Chib-style level-2 calculation | `mcmc_irt_hier1d(..., chib=.true.)`; result in `log_marginal` | `mcmcpack_irthier` |

`testpanelGroupBreak` and `testpanelSubjectBreak` are primarily R orchestration/model-comparison wrappers around fitted break/HMM models rather than independent low-level sampler engines. Their formula/indexing/presentation wrapper layer is not reproduced as a separate Fortran API; the underlying break and panel-HMM computational machinery is present.

`choicevar` is R formula-construction infrastructure and is intentionally not a separate Fortran routine. `read.Scythe`/`write.Scythe` are R/Scythe serialization helpers and are also omitted.

Plotting functions such as `plotChangepoint`, `plotHDPChangepoint`, `plotState`, `tomogplot`, and `dtomogplot` are intentionally omitted per the requested scope.

## Specialized-sampler implementation details

- **Hierarchical/dynamic EI:** retains the table-count EI likelihood and latent pair updates. The upstream adaptive/doubling slice mechanism is represented by a standard stepping-out/shrinkage slice transition targeting the same one-dimensional full conditionals.
- **Dynamic 1-D IRT:** uses Albert-Chib latent utilities, item Gaussian updates, scalar state-space forward-filter/backward-sampling for respondent trajectories, inverse-gamma innovation variances, and equality/inequality identification constraints.
- **Robust K-D IRT:** implements the four-parameter robust logistic item-response likelihood and contamination parameters. Coordinate random-walk Metropolis replaces the upstream coordinate slice transition while targeting the same posterior.
- **DP 2-D paired comparison:** implements a finite truncated stick-breaking DP for respondent-angle clusters, membership and concentration updates, angle moves, and constrained 2-D candidate-position updates.
- **SSVS quantile regression:** implements the asymmetric-Laplace latent-weight representation, inverse-Gaussian latent weights, exponential local shrinkage, Beta inclusion-probability update, and integrated Gaussian model-indicator probabilities. Covariance systems are recomputed directly rather than reproducing the C++ rank-one Cholesky/Givens optimization.
- **HMM panel FE/RE:** implements ordered hidden-state sampling and regime parameter updates. The RE state draw is conditional on currently sampled random effects rather than using the upstream collapsed marginal state probability; this is an alternative Gibbs blocking of the same joint model.
- **HDP HMMs:** use an explicit finite-`K` weak-limit hierarchical transition representation with CRT auxiliary counts and sticky self-transition mass. Poisson and negative-binomial emissions are updated with state-path sampling; negative-binomial dispersion uses the package parameterization and slice updates.
- **HDP HSMM:** uses explicit-duration segment sampling with zero-truncated negative-binomial durations and no self-transition. The current implementation conditions the terminal segment on ending at the observed final time; the upstream C++ performs a right-censored final-duration augmentation, so terminal-duration transition behavior is not bit-for-bit identical.
- **Hierarchical IRT PX/Chib:** includes the Liu-Wu-style parameter-expanded latent-data path and a Chib-style level-2 marginal-likelihood calculation exposed through the result type.

## Other deliberate numerical/kernel differences

Several MCMCpack C++ samplers use large hand-tuned auxiliary normal-mixture tables or Scythe-specific blocking. Where reproducing those tables would duplicate opaque implementation machinery, the Fortran code keeps the model/posterior target but uses a simpler valid transition:

- `MCMCnegbin`, `MCMCnegbinChange`, `HDPHMMnegbin`, and `HDPHSMMnegbin`: regression-coefficient updates use Gaussian random-walk Metropolis against the exact negative-binomial likelihood instead of the upstream auxiliary normal-mixture beta Gibbs step. Dispersion parameterization and slice updates follow the package model.
- `MCMCpoissonChange` and `HDPHMMpoisson`: regime beta updates use direct posterior Metropolis updates rather than the auxiliary normal-mixture construction.
- `MCMChregress`: uses algebraically equivalent conditional Gaussian Gibbs blocks rather than reproducing the exact Chib-Carlin block implementation.
- Marginal-likelihood/Chib calculations attached to several ordinary upstream R interfaces are not universally reproduced. Hierarchical IRT's specialized level-2 calculation is implemented in v0.2.0, and `bayes_factor` can combine supplied log marginal likelihoods elsewhere.

These differences alter Monte Carlo transition behavior and exact random streams, but the documented replacement updates are constructed for the same stated posterior target except for the explicitly noted HSMM terminal-duration boundary convention.

## Random-number streams

MCMCpack's C++ code supports Scythe Mersenne Twister and L'Ecuyer substreams. The Fortran translation uses the processor's standard `random_number` generator, with `set_seed` providing reproducible seeding on a given compiler/runtime. Draws are therefore not bit-for-bit reproducible against R MCMCpack and need not be identical across Fortran compilers.

## Dependency translations

The user-supplied translations are vendored and used directly:

- `coda-fortran`: sampler matrices can be converted with `as_coda_chain`.
- `quantreg-fortran`: `quantreg_start` supplies quantile-regression starting estimates.
- `mcmc-fortran`: generic Metropolis result/scaling types and `metrop` are re-exported through `mcmcpack_dependencies`.

Their original license files and source trees remain under `vendor/`.

## Validation

The release contains 22 test programs covering distributions/core utilities, ordinary regression samplers, hierarchical models, latent-variable models, factor/ordinal/mixed-factor models, multinomial logit, negative-binomial models, ordered changepoints, residual-break analysis, SVD regression, model utilities, EI/dynamic/robust IRT, DP paired comparison, SSVS quantile regression, panel HMMs, HDP-HMM/HDP-HSMM count models, and parameter-expanded hierarchical IRT with marginal-likelihood output.

A clean full-tree build of **58 Fortran source units** (package plus all three vendored dependencies) was performed with GNU Fortran 14.2.0 using standard free-form Fortran 2018 source-line limits:

```text
-std=f2018 -fcheck=all
```

All **22/22** test programs pass under runtime checking. No `-ffree-line-length-none` option is needed.

An FPM executable was not installed in the build environment. `fpm.toml` and all local path-dependency declarations are included, so the intended workflow is `fpm build` / `fpm test`; the equivalent compiler-level build and test run was performed directly.
