# Upstream API map

Fortran identifiers are case-insensitive. R methods that depended on S4
dispatch are represented by typed generic interfaces or explicit procedures.

| Upstream API | Fortran API | Status |
|---|---|---|
| `AICc` | `aicc` | Implemented |
| `V2b`, `V2p`, `V3b`, `V3p` | `v2b`, `v2p`, `v3b`, `v3p` | Implemented |
| `Vdegenerate`, `Vlinear`, `Vsymmetric` | `vdegenerate`, `vlinear`, `vsymmetric` | Implemented |
| `vtrans`, `vgradient`, `vinverse`, `vdownprob` | Same names | Implemented |
| `stochinverse`, `pcoincide` | Same names | Implemented |
| `gauss`, `gauss0` constructors | `margin('gauss',...)`, `margin('gauss0',...)` | Implemented |
| `laplace`, `laplace0` | `margin('laplace',...)`, `margin('laplace0',...)` | Implemented |
| skew Laplace | `margin('slaplace',...)` | Implemented |
| double/skew-double Weibull | `margin('doubleweibull',...)`, `margin('sdoubleweibull',...)` | Implemented |
| Student/centered/skew Student | `margin('st',...)`, `margin('st0',...)`, `margin('sst',...)` | Implemented |
| all exported `d/p/q/r` margin functions | Same names | Implemented |
| `dmarg`, `pmarg`, `qmarg` | Same names | Implemented |
| `edf`, `fitEDF`, `pedf`, `predict_empirical` | `edf`, `fit_edf`, `pedf`, `predict_empirical` | Implemented with Gaussian KDE |
| `armacopula` | `armacopula` / `arma_copula` | Implemented |
| `sarmacopula` | `sarmacopula` / `sarma_copula` | Implemented |
| `sarma2arma`, `expand_ar`, `expand_ma` | Same names | Implemented |
| `non_stat`, `non_invert` | Same names | Implemented |
| `starmaStateSpace` | internal `state_matrices` | Implemented internally |
| `kfilter` | `kfilter` | Implemented |
| `sigmastarma` | `sigmastarma` | Implemented |
| ARMA/SARMA `sim`, `predict`, `resid`, `fit` | generic `sim`, explicit prediction routines, `resid_arma_copula`, generic `fit` | Implemented |
| `armacopula_objective`, `sarmacopula_objective` | `arma_objective` after `sarma2arma` | Implemented |
| `acf2pacf`, `pacf2acf`, `pacf2ar` | Same names | Implemented |
| `kpacf_arma`, `kpacf_sarma4`, `kpacf_sarma12` | Same names | Implemented |
| `kpacf_arfima`, `kpacf_fbn` | Same names | Implemented by deterministic truncation |
| `glag`, `kendall` | `glag`, generic `kendall` | Implemented |
| `dvinecopula` | `dvinecopula` / `dvine_copula` | Implemented |
| `dvinecopula2` | `dvinecopula2` / `dvine2_copula` | Implemented |
| `dvinecopula3` | `dvinecopula3` / `dvine3_copula` | Implemented |
| `mklist_dvine`, `mklist_dvine2`, `mklist_dvine3` | Same names | Implemented |
| `Rblatt`, `IRblatt`, `Rblattdens` | `rblatt`, `irblatt`, `rblattdens` | Implemented |
| `simdvine` | `simdvine` | Implemented |
| D-vine objectives | `dvine_objective`, `dvine_loglik` | Implemented |
| D-vine `fit`, `predict`, `resid`, `kendall` | generic/typed interfaces | Implemented; fitting is sequential |
| `arma2dvine`, `armafit2dvine`, `sarma2dvine` | Same names | Implemented |
| `ktau_to_par` | `kendall_to_parameter` | Implemented |
| `swncopula` | `swncopula` | Implemented |
| `vtscopula` | `vtscopula` / `vtscopula_spec` | Implemented |
| `setwcopula`, `vtparlist` | Same names | Implemented |
| `vtscopula_objective`, `wobjective` | Same names | Implemented |
| `profilefulcrum` | `profilefulcrum` | Implemented numerically |
| `dcondvtarma`, `pcondvtarma`, `qcondvtarma` | Same names | Implemented |
| `tscm` | `tscm` / `tscm_spec` | Implemented |
| `fitSTEPS` | `fit_tscm_steps` / generic `fit` | Implemented |
| `fitFULLa`, `fitFULLb` | `fit_full` | Stepwise equivalent, not joint optimization |
| full-model `sim`, `predict`, `quantile`, `logLik` | `sim`, `predict_tscm_*`, `tscm_loglik` | Implemented |
| `safe_ses` | `safe_ses` | Implemented |
| `strank` | `strank` | Implemented |
| `setoptions` | Not applicable | R list-management utility omitted |
| plotting and S4 show/coerce methods | Not applicable | Omitted |
