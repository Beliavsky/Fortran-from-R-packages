# RobStatTM to Fortran API map

R names containing dots are normalized to valid Fortran identifiers. Compatibility wrappers omit dots and usually omit underscores; idiomatic implementations use snake case.

| Upstream R API | Preferred Fortran API | Compatibility API | Notes |
|---|---|---|---|
| `rho` | `rho_value` | `rho` | Scalar elemental loss. |
| `rhoprime` | `rho_prime` | `rhoprime` | First derivative. |
| `rhoprime2` | `rho_second` | `rhoprime2` | Second derivative. |
| `scaleM` | `scale_m` | `scalem` | M-scale fixed-point iteration. |
| `INVTR2` | `inverse_robust_r_squared` | `invtr2` | Inverts robust R-squared mapping. |
| `huber` | `tuning_huber` | `huber` | Gaussian-efficiency tuning. |
| `bisquare` | `tuning_bisquare` | `bisquare` | Gaussian-efficiency tuning. |
| `opt`, `optv0` | `tuning_opt` | `opt`, `optv0` | Scalar tuning interface. |
| `mopt`, `moptv0` | `tuning_mopt` | `mopt`, `moptv0` | Scalar tuning interface. |
| `locScaleM`, `MLocDis` | `loc_scale_m` | `mlocdis` | Robust location, scale, and standard error. |
| `MMPY` | `mm_py_fit` | `mmpy` | MM regression with robust S start. |
| `SMPY` | `sm_py_fit` | `smpy` | Numeric S/MM path without R factor metadata. |
| `lmrobM` | `lmrob_m_fit` | `lmrobm` | M regression result. |
| `lmrobdetMM` | `lmrobdet_mm` | `lmrobdetmm` | Matrix interface. |
| `DCML`, `lmrobdetDCML` | `dcml_fit`, `lmrobdet_dcml` | `lmrobdetdcml` | Direct and end-to-end forms. |
| `cov.dcml` | `dcml_covariance` | -- | DCML covariance kernel. |
| `lmrobdetMM.RFPE` | `robust_rfpe`, `lmrobdet_mm_rfpe` | `lmrobdetmm_rfpe` | Can return both RFPE terms. |
| `step.lmrobdetMM` | `stepwise_rfpe` | -- | Numeric column-level forward/backward/both selection. |
| `lmrobdetLinTest`, `rob.linear.test` | `robust_linear_test`, `lmrobdet_lin_test` | `lmrobdetlintest`, `roblineartest` | Nested robust models. |
| `BYlogreg`, `logregBY` | `logreg_by` | `bylogreg` | BY logistic estimator. |
| `WBYlogreg`, `logregWBY` | `logreg_wby` | `wbylogreg` | Leverage-screened BY estimator. |
| `WMLlogreg`, `logregWML` | `logreg_wml` | `wmllogreg` | Leverage-screened ML estimator. |
| `fastmve` | `fast_mve` | `fastmve` | rrcov-backed MVE numerical implementation. |
| `initPP`, `KurtSDNew` | `init_pp`, `kurt_sd_new` | `initpp` | Projection outlyingness initializer. |
| `covClassic` | `cov_classic` | `covclassic` | Classical covariance and distances. |
| `covRob`, `Multirobu` | `cov_rob` | `covrob`, `multirobu` | Automatic robust scatter dispatcher. |
| `covRobMM`, `MMultiSHR` | `cov_rob_mm` | `covrobmm`, `mmultishr` | MM-SHR scatter. |
| `covRobRocke`, `RockeMulti` | `cov_rob_rocke` | `covrobrocke`, `rockemulti` | Rocke scatter. |
| `SMPCA`, `pcaRobS` | `pca_rob_s`, `sm_pca` | `smpca` | Residual M-scale PCA. |
| `prcompRob` | `prcomp_rob` | `prcomprob` | Fortran result type analogous to `prcomp`. |

Result types are defined in `robstattm_types`: `regression_result`, `logistic_result`, `covariance_result`, `projection_result`, `location_scale_result`, `pca_result`, `linear_test_result`, and `model_selection_result`.
