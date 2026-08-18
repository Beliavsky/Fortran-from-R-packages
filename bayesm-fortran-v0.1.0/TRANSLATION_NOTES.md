# Translation notes

## Scope

This release translates the computational exports of `bayesm` 3.1-7 to
modern Fortran.  R plotting and S3 summary methods are intentionally omitted.
R lists are represented by derived types and allocatable arrays.

## Export mapping

| R export | Fortran |
|---|---|
| `breg` | `breg` |
| `createX` | `create_x` |
| `eMixMargDen` | `e_mix_marg_den` |
| `mixDen`, `mixDenBi`, `momMix` | `mix_den`, `mix_den_bi`, `mom_mix` |
| `llmnl`, `mnlHess` | `llmnl`, `mnl_hess` |
| `llmnp`, `mnpProb`, `ghkvec` | `llmnp`, `mnp_prob`, `ghkvec` |
| `llnhlogit`, `simnhlogit` | `llnhlogit`, `simnhlogit` |
| `lndIChisq`, `lndIWishart`, `lndMvn`, `lndMvst` | `lnd_ichisq`, `lnd_iwishart`, `lnd_mvn`, `lnd_mvst` |
| `nmat`, `numEff`, `condMom`, `logMargDenNR`, `cgetC` | `nmat`, `num_eff`, `cond_mom`, `log_marg_den_nr`, `cget_c` |
| `rdirichlet`, `rwishart`, `rmvst`, `rtrun` | `rdirichlet`, `rwishart`, `rmvst`, `rtrun` |
| `runireg`, `runiregGibbs`, `rmultireg`, `rsurGibbs` | `runireg`, `runireg_gibbs`, `rmultireg`, `rsur_gibbs` |
| `rbprobitGibbs`, `rordprobitGibbs` | `rbprobit_gibbs`, `rordprobit_gibbs` |
| `rmnpGibbs`, `rmvpGibbs`, `rbiNormGibbs` | `rmnp_gibbs`, `rmvp_gibbs`, `rbinorm_gibbs` |
| `rmnlIndepMetrop` | `rmnl_indep_metrop` |
| `rnegbinRw`, `rhierNegbinRw` | `rnegbin_rw`, `rhier_negbin_rw` |
| `rmixture`, `rmixGibbs`, `rnmixGibbs`, `clusterMix` | `rmixture`, `rmix_gibbs`, `rnmix_gibbs`, `cluster_mix` |
| `rhierLinearModel`, `rhierLinearMixture` | `rhier_linear_model`, `rhier_linear_mixture` |
| `rhierMnlRwMixture`, `rhierBinLogit` | `rhier_mnl_rw_mixture`, `rhier_bin_logit` |
| `rDPGibbs`, `rhierMnlDP` | `rdp_gibbs`, `rhier_mnl_dp` |
| `rivGibbs`, `rivDP` | `riv_gibbs`, `riv_dp` |
| `rscaleUsage` | `rscale_usage` |
| `rbayesBLP` | `rbayes_blp` |

## Numerical implementation

The package is self-contained.  Rcpp/Armadillo operations are replaced by
native Fortran Cholesky/SPD solves, matrix operations, RNGs, and probability
functions.  MCMC outputs are stored in explicit derived types instead of R
lists.

The regression, SUR, finite-normal-mixture, probit, MNL, negative-binomial,
hierarchical-linear, hierarchical finite-mixture, IV-Gibbs, scale-usage, and
BLP blocks retain the corresponding Bayesian conditional/MH structure of the
upstream algorithms.

### Dirichlet-process routines

`rdp_gibbs` is a native normal Dirichlet-process mixture sampler using CRP
allocation updates, normal/inverse-Wishart component updates, and an
Escobar-West concentration update.  The upstream Rcpp code parameterizes its
base-measure/hyperparameter grid somewhat differently, so this is a
posterior-model-equivalent DP implementation rather than an
instruction-for-instruction port.

`rhier_mnl_dp` uses the standard finite weak-limit representation
`Dirichlet(alpha/K, ..., alpha/K)` and approaches the DP prior as `K` grows.
It reuses the directly translated hierarchical finite-mixture MNL transition.

`riv_dp` currently preserves the directly translated correlated-error IV Gibbs
block and the DP concentration parameter/output plumbing.  Observation-level
DP residual density estimation is available separately through `rdp_gibbs`;
this v0.1 interface does not reproduce every latent residual-mixture update in
upstream `rivDP_rcpp_loop.cpp`.

### Scale-usage routine

The ordinal latent-utility, threshold, covariance, mean, and GHK-based
classification transitions are translated.  The Lambda sufficient-statistic
update uses the sample covariance representation rather than the exact
upstream helper's internal moment bookkeeping.

## Deliberately omitted R infrastructure

- `plot.bayesm.*`
- `summary.bayesm.*`
- R printing/classes and nested-list conventions
- R progress printing and console flushing

The original R and C++ sources are retained under `upstream/` for provenance,
but they are not compiled or linked.
