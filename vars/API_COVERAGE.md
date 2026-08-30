# API coverage

Legend: **covered** = numerical behavior has a native Fortran entry point;
**array-level** = computation is covered but R object extraction/assembly is not;
**external** = upstream delegates the computation to another package;
**omitted interface** = presentation, formula, S3, or callback machinery.

| R API / area | Status | Fortran API / disposition |
| --- | --- | --- |
| `VAR` | covered | `fit_var` |
| `VARselect` | covered | `var_select` |
| `restrict` manual | covered | `restrict_var_manual` |
| `restrict` SER | covered | `restrict_var_ser` |
| `A`, `Acoef`, `B`, `Bcoef` | covered | typed `var_model%a`, `%coef`, `structural_impact` |
| `Phi` | covered | `phi_from_a`, `structural_impulse_response` uses the same recursion |
| `Psi` | covered | `psi_from_a_sigma` |
| `roots` | covered | `var_roots` |
| `predict.varest` | covered | `forecast_var`; future exogenous/non-lag regressors can be supplied explicitly |
| forecast covariance | covered | `forecast_covariance` |
| `irf.varest` point IRF | covered | `impulse_response` |
| `irf.varest` residual bootstrap | covered | `bootstrap_irf_indices`; caller supplies resampling indices |
| `irf.svarest` / `irf.svecest` point IRF | covered | `structural_impulse_response` |
| structural/SVEC bootstrap re-estimation | partial | structural estimators and bootstrap path pieces are available separately; R object/RNG orchestration is not cloned |
| `fevd.varest` | covered | `fevd_var` |
| structural/SVEC FEVD | covered | `structural_fevd` |
| `BQ` | covered | `bq_identification` |
| `arch.test` | covered | `arch_test_univariate`, `arch_test_multivariate` |
| `normality.test` | covered | `jarque_bera_univariate`, `jarque_bera_multivariate` |
| `serial.test` PT/PT.adjusted/BG/ES | covered | `portmanteau_tests`, `bg_serial_tests` |
| `causality` Granger | covered (unrestricted VAR) | `granger_causality`; coefficient-restricted R-model overlays are not reproduced |
| `causality` instantaneous | covered | `instantaneous_causality` |
| optional `vcov.`/sandwich callback in causality | omitted interface | callback into R model/sandwich infrastructure; no duplicated dependency code |
| `SVAR(..., estmethod="scoring")` | covered | `svar_fit_scoring` |
| `SVAR(..., estmethod="direct")` | computational objective covered | `svar_negloglik`; generic optimizer is intentionally caller-selected |
| SVAR LR calculation | covered | returned in `svar_result` when identified/overidentified |
| `SVEC` long-run multiplier | covered | `svec_long_run_matrix` |
| `SVEC` scoring | array-level | `svec_fit_scoring` accepts alpha, beta, gamma, covariance and restriction masks |
| `vec2var` lag conversion | array-level | `vec2var_coefficients` for `transitory` and `longrun` specifications |
| extraction from `urca::ca.jo` | external/object interface | use numerical arrays from the sibling translated `urca` package; no `urca` source is vendored |
| `logLik.varest` | covered | `var_loglik` |
| structural log-likelihood objective | covered | `svar_negloglik` |
| `stability` | external | upstream calls `strucchange::efp`; use sibling translated `strucchange` package |
| `coef`, `fitted`, `residuals` methods | covered as data | fields of `var_model`; S3 wrappers omitted |
| `toMlm`, `coeftest`, `bread`, `estfun`, `vcovHC` integration | omitted interface | R model/sandwich interoperability rather than core VAR numerics |
| print/summary/plot/fanchart | omitted interface | presentation-only R code |

## Deterministic test coverage

`test/test_vars.f90` checks a fixed bivariate VAR against independent numerical
references for coefficient estimates, lag selection, MA recursion, roots,
forecasts, covariance, BQ identification, FEVD, univariate and multivariate
normality, ARCH, Portmanteau, BG and Edgerton-Shukur diagnostics, Granger and
instantaneous causality, identity-resample bootstrap reconstruction, both `vec2var`
specifications, the SVEC long-run multiplier and scoring estimator, and SVAR
scoring/likelihood.
