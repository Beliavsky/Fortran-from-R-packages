# Numerical API map

| Upstream concept | Fortran API | Status |
|---|---|---|
| GEV d/p/q/r | `dgev`, `pgev`, `qgev`, `rgev` | translated |
| GPD d/p/q/r | `dgp`, `pgp`, `qgp`, `rgp` | translated |
| Extended GPD d/p/q/r | `degp`, `pegp`, `qegp`, `regp` | all 7 transforms |
| GEV/GPD likelihood | `gev_ll`, `gpd_ll` | translated |
| GEV/GPD score/information | `gev_score`, `gpd_score`, `gev_infomat`, `gpd_infomat` | translated |
| GEV/GPD MLE | `gev_fit`, `gpd_fit` | translated |
| EGP fit/return level | `egp_fit`, `egp_retlev` | translated |
| GEV return level | `gev_retlev`, `gev_nyr` | translated |
| N-block parameterizations | `gpd_n_mean`, `gpd_n_quant`, `gev_n_mean`, `gev_n_quant` | translated |
| Expected-shortfall GPD likelihood | `gpde_ll`, `gpde_score`, `gpde_infomat` | translated |
| Return-level GPD likelihood | `gpdr_ll`, `gpdr_score`, `gpdr_infomat` | translated |
| Return-level GEV likelihood | `gevr_ll`, `gevr_score`, `gevr_infomat` | translated |
| N-block GPD/GEV likelihood | `gpdn_*`, `gevn_*` | translated |
| r-largest order statistics | `rrlarg`, `rlarg_ll`, `rlarg_score`, `rlarg_infomat` | translated |
| Poisson-process likelihood | `pp_ll`, `pp_score`, `pp_infomat` | translated |
| Tail-index estimators | `shape_hill`, `shape_pickands`, `shape_moment`, `shape_vries`, `shape_genjack`, `shape_osz`, `shape_genquant` | translated |
| Pickands-Xu | `pickands_xu` | translated |
| Second-order rho estimators | `rho_dk`, `rho_fagh`, `rho_ghp`, `rho_gbw` | translated |
| Weissman quantile/CI | `qweissman`, `qweissman_ci` | translated |
| PWM/L-moments | `pwm`, `lmoments`, `gpd_lmom` | translated |
| Mean residual life | `mrl_profile` | translated |
| Semiparametric margins | `spunif_vector`, `spunif_matrix` | translated, explicit tail params |
| Tail dependence | `taildep_empirical`, `taildep_hill` | translated |
| Extremogram | `xacf_extremogram` | translated |
| Extremal coefficients | `extcoef_fmado`, `extcoef_smith`, `extremo_pairwise` | translated |
| Euclidean likelihood | `euclidean_weights` | translated |
| Empirical likelihood | `empirical_likelihood` | translated |
| Empirical Pickands | `pickands_empirical` | translated |
| Dirichlet log density | `ldirfn` | translated |
| Anisotropic distance | `distg`, `dgeoaniso` | translated |
| Variogram/correlation | `power_vario`, `powerexp_cor`, `schlather_vario` | translated |
| Huesler-Reiss covariance transform | `lambda2cov` | translated |
| Dirichlet RNG | `rdir` | translated |
| MV normal/t RNG | `rmnorm`, `mvrt_sample`, `dmvnorm` | translated |
| Spectral RNGs | `rlogspec`, `rneglogspec`, `rbilogspec`, `rdirspec`, `rhrspec`, `rdirmixspec` | translated |
| Generic spectral RNG | `rmev_spectral` | translated for listed models |
| Exact max-stable Algorithm 1 | `rmev` | translated for listed models |
| GPD/Pareto transforms | `gpd_to_pareto`, `jac_gpd_pareto`, `ordexp_to_gev` | translated |
| Lower/trimmed Hill | `shape_lthill`, `shape_lthill_path`, `shape_trimhill` | translated |
| BAB threshold selection | `bab_fcst`, `thselect_bab` | translated; uses vendored `expint` |
| GPD Cox-Snell/Firth bias correction | `gpd_bias`, `gpd_fscore`, `gpd_bcor` | translated; uses vendored `nleqslv` |
| Stein weighted GPD | `stein_weights`, `stein_gp_lik`, `fit_wgpd` | translated |
| Exponential-regression tail index | `shape_erm` (`bdgm`, `fh`) | translated |
| Krupskii-Joe tail dependence | `kjtail`, `kjtail_uniform` | translated |
| GPD/GEV profile likelihood | `gpd_profile`, `gev_profile` | translated for common fixed/reparameterized paths |
| MV normal/t upper probabilities | `mvn_upper_prob_qmc`, `mvt_upper_prob_qmc` | native deterministic QMC |
| Multivariate exponent measures | `expme_logistic`, `expme_neglog`, `expme_br`, `expme_br_wt`, `expme_hr`, `expme_xstud` | translated |
| GEV Cox-Snell/Firth bias correction | `gev_bias`, `gev_fscore`, `gev_bcor` | translated; deterministic cumulant quadrature + nleqslv |
| BAB Monte Carlo test/envelope | `thselect_bab(..., test=.true.)` | translated |
| GPD/GEV TEM profile correction | `gpd_tem_profile`, `gev_tem_profile` | translated for principal original parameters |
| Multivariate GP likelihood | `mgp_ll_log`, `mgp_ll_neglog`, `mgp_ll_br`, `mgp_ll_xstud` | translated |
| Censored multivariate GP likelihood | `mgp_cll_log`, `mgp_cll_neglog`, `mgp_cll_br`, `mgp_cll_xstud` | translated; native Gaussian/t QMC |
| Extremal-t spectral RNG | `rexstudspec` | translated |
| Brown-Resnick spectral RNG | `rbrspec` | translated (covariance parameterization) |
| Extended max-stable Algorithm 1 | `rmev` | now includes extremal-t and Brown-Resnick |
| R-Pareto simulation | `rparp` | translated by direct rejection for common risk functionals |
| Generalized R-Pareto simulation | `rgparp` | translated by direct rejection for common risk functionals |
