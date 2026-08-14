# R-to-Fortran API mapping

The Fortran interface is intentionally array/type based rather than a copy of
R's S3 object and formula conventions.

| Upstream R family | Fortran API | v0.1.2 status |
|---|---|---|
| `em*`, `me*` | `fit_model` | Implemented for E/V and all 14 multivariate models |
| `mstep*` | `mstep_model` | Implemented |
| `estep*`, `cdens*` | `mixture_posterior`, `component_log_density`, `mixture_log_density` | Implemented |
| `mvnX`, `mvnXII`, `mvnXXI`, `mvnXXX` | one-component E/EII/EEI/VVV fits and density helpers | Computational equivalent |
| `sim*` | `simulate_mixture`, `simulate_fit` | Implemented |
| `nVarParams`, `nMclustParams` | `n_var_params`, `n_mclust_params` | Implemented |
| `bic`, `pickBIC`, `mclustBIC` | `bic_value`, `pick_bic`, `mclust_bic` | Implemented |
| `Mclust` | `mclust_select` | Implemented for ordinary Gaussian mixtures |
| `icl`, `mclustICL` | `icl_value` / selected fit `%icl` | Implemented |
| `hc`, `hcE`, `hcV`, `hcEII`, `hcVII`, `hcEEE`, `hcVVV` | `hc_fit`, `hc_responsibilities`, `hclass` | Implemented |
| `map` | `map_z` | Implemented |
| `unmap` | `unmap_classes` | Implemented |
| `matchCluster` | `match_clusters` | Implemented with maximum-overlap assignment |
| `majorityVote` | `majority_vote` | Implemented |
| `adjustedRandIndex` | `adjusted_rand_index` | Implemented |
| `BrierScore` | `brier_score` | Implemented |
| `logsumexp`, `softmax` | same names | Implemented |
| `covw` | `covariance_weighted` | Implemented |
| `dmvnorm` | `dmvnorm` | Implemented |
| `hdrlevels` | `hdr_levels` | Implemented |
| `count` | `count_values` | Implemented |
| `densityMclust`, `dens` | density helpers on `mclust_fit` | Implemented |
| `cdfMclust`, `quantileMclust` | `cdf_mclust_1d`, `quantile_mclust_1d` | 1-D implementation |
| `MclustDA` | `fit_mclust_da`, `predict_mclust_da` | Implemented |
| `MclustSSC` | `fit_model_ssc`, `mclust_ssc_select` | Implemented |
| `MclustDR` | `fit_mclust_dr`, `project_mclust_dr` | Implemented |
| `crimcoords` | `fit_crimcoords` | Implemented |
| `imputeData` | `impute_data` | Implemented |
| `me.weighted` | `fit_model_weighted` | Implemented |
| `clustCombi` | `clust_combi`, `cluster_combination`, `apply_combination` | Implemented core entropy merging |
| `mclustBootstrapLRT` | `bootstrap_lrt` | Implemented parametric bootstrap core |
| `MclustBootstrap` | `mclust_parameter_bootstrap` | Implemented parametric parameter bootstrap core |
| `simulate.Mclust`, `simulate.densityMclust` | `simulate_fit` | Implemented |
| `randomOrthogonalMatrix` | `random_orthogonal_matrix` | Implemented |
| `priorControl`, `defaultPrior`, prior `*p` paths | -- | Not in v0.1.2 |
| background/noise-component BIC fits | -- | Not in v0.1.2 |
| `gmmhd*` | -- | Not in v0.1.2 |
| `MclustDRsubsel*` | -- | Not in v0.1.2 |
| `cvMclustDA` | -- | Full wrapper not in v0.1.2 |
| `imputePairs` | -- | Not in v0.1.2 |
| `hcRandomPairs`, `randomPairs`, `dupPartition` | -- | Not in v0.1.2 |
| plotting/print/summary/S3 conversion functions | -- | Intentionally omitted interface/UI code |

The umbrella module is `mclust`; applications generally need only `use mclust`.
