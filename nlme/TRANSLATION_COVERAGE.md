# Translation coverage

The original namespace contains more than one hundred exported names, most of
which are constructors, S3 methods, formula/data utilities, or display methods.
This project translates the numerical model layer and provides typed Fortran
replacements for the principal computational exports.

| Original area | Fortran representation | Coverage |
|---|---|---|
| `gls`, `glsEstimate`, `glsApVar` | `fit_gls`, result covariance and likelihood fields | Implemented |
| `lme`, `fixed.effects`, `random.effects`, `getVarCov`, `VarCorr` | `fit_lme`, `lme_result` | Implemented |
| `gnls` | `fit_gnls` | Implemented with fixed covariance parameters |
| `nlme` | `fit_nlme` | Implemented by first-order linearization |
| `corAR1`, `corCAR1`, `corARMA`, `corCompSymm` | `correlation_spec`, `correlation_matrix` | Implemented |
| `corExp`, `corGaus`, `corLin`, `corRatio`, `corSpher` | Spatial correlation kinds | Implemented |
| `corSymm`, `corNatural` | `COR_UNSTRUCTURED` normalized Cholesky | Adapted |
| `corFactor`, `corMatrix`, `logDet` | Matrix construction and SPD helpers | Implemented |
| `varFixed`, `varIdent`, `varPower`, `varExp` | `variance_spec` | Implemented |
| `varConstPower`, `varConstProp` | `variance_spec` | Implemented |
| `varComb` | Precombine covariates/scales or construct covariance explicitly | Not a separate container |
| `pdIdent`, `pdDiag`, `pdCompSymm`, `pdLogChol` | `pd_spec` | Implemented |
| `pdSymm`, `pdNatural` | `PD_LOG_CHOL` | Adapted |
| `pdBlocked`, `reStruct` | Explicit `Z` construction | Not a separate container |
| `ACF`, `Variogram`, `pooledSD` | Diagnostic procedures | Implemented |
| `simulate.lme` | `simulate_lme` | Implemented |
| `lmList`, `nlsList` | `fit_lm_list`, `fit_nls_list` | Implemented |
| `gsummary`, `isBalanced` | `group_summary`, `is_balanced` | Implemented |
| `fdHess` | `finite_difference_hessian` | Implemented |
| Formula methods, S3 accessors/replacements, grouped-data classes | Explicit arrays and result fields | Replaced by typed API |
| Plotting, `augPred`, `comparePred`, lattice methods | None | Omitted |
| Package datasets and scripts | None | Omitted from the compiled/source translation |
