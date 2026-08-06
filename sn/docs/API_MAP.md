# API map

## Univariate distributions

| R API | Fortran API | Module |
|---|---|---|
| `dsn`, `psn`, `qsn`, `rsn` | same lowercase names | `sn_univariate` |
| `dst`, `pst`, `qst`, `rst` | same lowercase names | `sn_univariate` |
| `dsc`, `psc`, `qsc`, `rsc` | same lowercase names | `sn_univariate` |
| `T.Owen` | `owen_t` | `sn_univariate` |
| `zeta` | `zeta` | `sn_univariate` |
| `sn.cumulants` | `sn_cumulants` | `sn_univariate` |
| `st.cumulants` | `st_cumulants` | `sn_univariate` |
| `modeSECdistr` | `mode_sn`, `mode_st`, `mode_sc` | `sn_univariate` |
| `fournum` | `fournum` | `sn_univariate` |
| `dp2cp`, `cp2dp` for SN | `dp_to_cp_sn`, `cp_to_dp_sn` | `sn_univariate` |
| `dp2op`, `op2dp` | `dp_to_op_uv`, `op_to_dp_uv` | `sn_parameters` |

## Multivariate distributions

| R API | Fortran API | Module |
|---|---|---|
| `dmsn`, `pmsn`, `rmsn` | same lowercase names | `sn_multivariate` |
| `dmst`, `pmst`, `rmst` | same lowercase names | `sn_multivariate` |
| `dmsc`, `pmsc`, `rmsc` | same lowercase names | `sn_multivariate` |
| `marginalSECdistr` | `marginal_sn` | `sn_multivariate` |
| `affineTransSECdistr` | `affine_transform_sn` | `sn_multivariate` |
| `conditionalSECdistr` | `conditional_sn` | `sn_multivariate` |
| multivariate `dp2op`, `op2dp` | `dp_to_op_mv`, `op_to_dp_mv` | `sn_parameters` |
| `msn.mle`, `msn.mple` | `msn_mle`, `msn_mple` | `sn_fit` |
| `mst.mple` | `mst_mple` | `sn_fit` |

The multivariate fit routines are deterministic marginal-composite estimators,
not exact replicas of the upstream joint optimizer.

## SUN

| R API | Fortran API | Module |
|---|---|---|
| `dsun`, `psun`, `rsun` | same lowercase names | `sn_sun` |
| `sunMean`, `sunVcov` | `sun_moments` fields | `sn_sun` |
| `marginalSUNdistr` | `marginal_sun` | `sn_sun` |
| `affineTransSUNdistr` | `affine_transform_sun` | `sn_sun` |
| `conditionalSUNdistr` | `conditional_sun_equal`, `conditional_sun_greater` | `sn_sun` |
| `joinSUNdistr` | `join_sun` | `sn_sun` |
| `convolutionSUNdistr` | `convolution_sun` | `sn_sun` |
| `convertSN2SUNdistr` | `sn_to_sun` | `sn_sun` |
| `convertCSN2SUNpar` | `csn_to_sun` | `sn_sun` |

## Fitting and utilities

| R API | Fortran API | Module |
|---|---|---|
| `selm`, `selm.fit` | `selm_fit` | `sn_fit` |
| `sn.mple`, `st.mple` | `sn_mple`, `st_mple` | `sn_fit` |
| `predict.selm` | `predict_selm` | `sn_fit` |
| `fitdistr.grouped` | `fit_grouped` | `sn_fit` |
| `Qpenalty` | `q_penalty` | `sn_misc` |
| `dSymmModulated` | `d_symm_modulated` | `sn_misc` |
| `rSymmModulated` | `r_symm_modulated` | `sn_misc` |
| `pprodn2`, `pprodt2`, `qprodt2` | same lowercase names | `sn_misc` |
| `galton_moors2alpha_nu` | `galton_moors_to_alpha_nu` | `sn_misc` |
| `vech`, `vech2mat` | `vech`, `vech_to_matrix` | `sn_matrix` |
| `duplicationMatrix` | `duplication_matrix` | `sn_matrix` |
| `tr`, `blockDiag` | `trace`, `block_diag` | `sn_matrix` |

## Intentionally omitted R interfaces

Plotting methods, S3/S4 constructors and display methods, formulas, model frames,
profile-likelihood plotting, R data-frame handling, and namespace hooks are not
computational Fortran APIs. Exact `MPpenalty`, the `sn.info*`/`st.info*` symbolic
information objects, and `quantreg`-based preliminary fitting are not ported.
