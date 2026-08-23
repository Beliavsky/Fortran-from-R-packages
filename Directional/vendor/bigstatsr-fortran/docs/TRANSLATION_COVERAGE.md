# Translation coverage

Upstream package: `bigstatsr` 1.6.2.

## Direct numerical equivalents

| R export / area | Fortran API | Status |
|---|---|---|
| `AUC` | `auc` | translated |
| `AUCBoot` | `auc_bootstrap` | translated |
| `FBM` | `fbm_real`, `create_fbm`, `attach_fbm` | translated for real64 |
| `FBM.code256` | `fbm_code256` | translated storage/count core |
| `as_FBM` | `create_fbm` + `fbm_from_array` | translated equivalent |
| `big_colstats` | `big_colstats` | translated |
| `big_counts` | `big_counts_rows`, `big_counts_cols` | translated |
| `big_prodVec` | `big_prod_vec` | translated |
| `big_cprodVec` | `big_cprod_vec` | translated |
| `big_prodMat` | `big_prod_mat` | translated |
| `big_cprodMat` | `big_cprod_mat` | translated |
| `big_crossprodSelf` | `big_crossprod_self` | translated |
| `big_tcrossprodSelf` | `big_tcrossprod_self` | translated |
| `big_cor` | `big_cor` | translated |
| `big_scale` | `big_scale` | translated |
| `big_transpose` | `fbm_transpose` | translated |
| `big_copy` | `fbm_copy` / array conversion | core translated |
| `big_increment` | `fbm_increment` | scalar/subset core translated |
| `big_univLinReg` | `big_univ_linreg` | translated |
| `big_univLogReg` | `big_univ_logreg` | translated; no R `glm` fallback |
| `big_spLinReg` | `elastic_net_gaussian_path` | numerical path translated |
| `big_spLogReg` | `elastic_net_logistic_path` | numerical path translated |
| `big_SVD` | `big_svd` | translated |
| `big_randomSVD` | `big_random_svd` | translated through ARPACK |
| SVD prediction | `svd_predict` | translated |
| sparse-model prediction | `predict_enet` | translated |
| `get_beta` | `get_beta` | translated |
| `pcor` | `pcor` | translated for numeric covariates |
| `block_size` | `block_size` | translated |
| `rows_along`, `cols_along` | same names | translated |
| biglasso summary kernel | `big_summaries` | translated |

## Functionality represented differently

- `add_code256`: the Fortran coded matrix is already an explicit type; callers
  supply the integer code mapping used for counts rather than wrapping an R FBM.
- `big_apply` / `big_parallelize`: numerical routines are block-oriented
  internally, but R's arbitrary-function/list-combine/cluster API is not a
  meaningful direct Fortran abstraction.
- `predict.big_sp`, `predict.big_SVD`: numerical prediction is available as
  `predict_enet` and `svd_predict`; S3 dispatch is not reproduced.

## R/runtime or presentation functionality intentionally omitted

- `big_attach`, `big_attachExtdata`, RDS metadata and reference classes
- `big_read` and `big_write` front ends that delegate to optional R packages
  `bigreadr` / `data.table`
- `covar_from_df` factor/data-frame encoding
- `nb_cores`, R cluster registration, `foreach`, `plus`
- `as_scaling_fun` closure construction
- `sub_bk`, `pasteLoc`, `without_downcast_warning`
- `asPlotlyText`, `plot_grid`, `theme_bigstatsr`
- all S3/S4 plotting, summary, extraction, replacement and printing methods

## Remaining worthwhile numerical parity targets

1. additional FBM storage encodings (`float`, `integer`, `unsigned short`)
2. the complete CMSA/cross-validation/grid-search wrapper used by
   `big_spLinReg` and `big_spLogReg`
3. OpenMP/block-parallel execution policy equivalent to the C++ kernels
4. a native portable delimited-file reader/writer if text I/O is desired
