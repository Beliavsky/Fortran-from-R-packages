# API map

| R function | Fortran API | Coverage |
|---|---|---|
| `wt.var` | `wt_var`, `weighted_variance` | Direct |
| `wt.moments` | `wt_moments`, `weighted_moments` | Direct |
| `wt.scale` | `wt_scale`, `weighted_scale` | Direct; metadata is returned in a type |
| `estimate.lambda` | `estimate_lambda` | Direct formula with smaller-Gram optimization |
| `estimate.lambda.var` | `estimate_lambda_var` | Direct |
| `var.shrink` | `var_shrink`, `variance_shrinkage` | Direct |
| `cor.shrink` | `cor_shrink`, `correlation_shrinkage` | Direct |
| `cov.shrink` | `cov_shrink`, `covariance_shrinkage` | Direct |
| `invcor.shrink` | `invcor_shrink`, `inverse_correlation_shrinkage` | Direct |
| `invcov.shrink` | `invcov_shrink`, `inverse_covariance_shrinkage` | Direct |
| `pcor.shrink` | `pcor_shrink`, `partial_correlation_shrinkage` | Direct; SPV returned in result type |
| `pvar.shrink` | `pvar_shrink`, `partial_variance_shrinkage` | Direct |
| `powcor.shrink` | `powcor_shrink`, `correlation_power_shrinkage` | Direct low-rank identity |
| `crossprod.powcor.shrink` | `crossprod_powcor_shrink` | Direct low-rank product |
| `cor2pcor` | `cor2pcor`, `correlation_to_partial` | Direct |
| `pcor2cor` | `pcor2cor`, `partial_to_correlation` | Direct |
| `mpower` | `mpower`, `matrix_power` | Direct symmetric eigendecomposition |
| `fast.svd` | `fast_svd` | Compact SVD using the smaller Gram matrix |
| `pseudoinverse` | `pseudoinverse` | Direct Moore-Penrose construction |
| `rank.condition` | `rank_condition` | Direct |
| `is.positive.definite` | `is_positive_definite` | Eigenvalue method |
| `make.positive.definite` | `make_positive_definite` | Higham eigenvalue correction |
| `rebuild.cov` | `rebuild_cov` | Direct |
| `decompose.cov` | `decompose_cov` | Direct |
| `rebuild.invcov` | `rebuild_invcov` | Direct |
| `decompose.invcov` | `decompose_invcov` | Direct |
| `sm2vec` | `sm2vec` | Direct lower-triangle ordering |
| `sm.index` | `sm_index` | Direct 1-based indices |
| `vec2sm` | `vec2sm` | Direct; missing diagonal represented by IEEE NaN |
