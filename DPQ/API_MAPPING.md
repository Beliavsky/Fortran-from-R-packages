# DPQ R-to-Fortran API mapping

Names are converted to conventional Fortran snake_case. Optional R
`lower.tail` and `log.p` arguments become optional logical arguments where
applicable.

| Upstream computational area / R names | Fortran module / representative API |
|---|---|
| `.D_*`, `.DT_*`, DPQ probability macros | `dpq_core`: `d_*`, `dt_*`, `prob_from_input`, `prob_output` |
| `log1mexp`, `log1pexp`, `logspace.add/sub` | `dpq_core`: `log1mexp`, `log1pexp`, `logspace_add/sub` |
| `pow*`, `expm1x*`, `p1l1*` | `dpq_core`: `dpq_pow`, `pow_di`, `pow1p`, `expm1x`, `expm1x_tser`, `p1l1*` |
| `logcf*`, `log1pmx`, `lgamma1p*`, `logr` | `dpq_core` |
| `lsum`, `lssum`, Chebyshev routines | `dpq_core` |
| `frexp`, `ldexp`, `modf`, `dpsifn` | `dpq_core`: `dpq_frexp`, `dpq_ldexp`, `dpq_modf`, `dpsifn_scalar` |
| `rexpm1`, `rlog1` | `dpq_core` |
| `bd0*`, `ebd0`, `stirlerr*`, `lgammacor` | `dpq_gamma_discrete` |
| `dpois_raw`, `dpois_simpl*`, `dgamma.R` | `dpq_gamma_discrete` |
| `dbinom_raw`, `dnbinomR`, `dnbinom.mu` | `dpq_gamma_discrete` |
| `ppoisD`, `ppoisErr`, `qpoisR`, `qbinomR`, `qnbinomR` | `dpq_gamma_discrete` |
| `algdiv`, `bpser`, `gam1*`, `gamln1*` | `dpq_gamma_discrete` |
| `qchisqAppr`, `qchisqKG`, `qchisqWH` | `dpq_gamma_discrete`: `qchisq_appr`, `qchisq_kg`, `qchisq_wh` |
| `qgammaAppr*`, `qgamma.R` | `dpq_gamma_discrete`: `qgamma_appr`, `qgamma_appr_kg`, `qgamma_appr_smallp`, `qgamma_r` |
| `pnormU_S53`, `pnormL_LD10`, `pnormAsymp` | `dpq_normal_beta` |
| `qnorm*` | `dpq_normal_beta`: `qnorm_r`, `qnorm_uappr*`, `qnorm_appr`, `qnorm_asymp`, `qnorm_cappr` |
| `qbetaAppr*`, `qbeta.R` | `dpq_normal_beta` |
| `pbetaAS_eq20/21`, `pbetaNorm2`, `pbetaRv1` | `dpq_normal_beta` |
| `pnbetaAppr2*`, `pnbetaAS310` | `dpq_normal_beta` |
| `lbeta_*`, `betaI`, `logQab_asy`, `Qab_terms` | `dpq_normal_beta` |
| Hypergeometric approximation family | `dpq_hyper`: `phyper_appr_as152`, `phyper_ibeta`, `phyper*_molenaar`, `phyper_peizer`, `phyper_bin*` |
| `pdhyper`, `phyperR`, `phyperR2` | `dpq_hyper`: `pdhyper`, `phyper_r`, `phyper_r2` |
| `dnchisqR`, `dnchisqBessel`, `dnoncentchisq` | `dpq_nchisq` |
| `pnchisq*`, `pnchi1sq`, `pnchi3sq` | `dpq_nchisq` |
| Patnaik/Pearson/Abdel-Aty/Sankaran/BolKuz/T93 | `dpq_nchisq`: corresponding snake_case functions |
| `qnchisq*`, `qchisqAppr.0-.3`, CF approximations | `dpq_nchisq` |
| `r_pois` | `dpq_nchisq:r_pois` |
| `b_chi*`, `lb_chi*`, `c_dt*`, `c_pt` | `dpq_t` |
| `pntR`, `pntR1`, `pntLrg`, `pntJW39*` | `dpq_t` |
| `pnt3150`, `pntP94`, `pntChShP94`, `pntVW13`, `pntGST23_*` | `dpq_t`: corresponding snake_case functions; see consolidation note |
| `dntJKBf*`, `dtWV` | `dpq_t` |
| `qtR*`, `qtAppr`, `qtNappr`, `qtU*`, `qntR*` | `dpq_t` |
| Wiener-germ `h*`, `g*`, `sW`, `qs`, `z*`, `pchisqW*` | `dpq_wiener` |
| `lgammaP11`, `dltgammaInc` | `dpq_toms1006` |

## Not translated as Fortran numerical APIs

- `pl2curves`, `plRpois`: plotting.
- `format01prec`: formatting/presentation.
- `all_mpfr`, `any_mpfr`, Rmpfr dispatch: R type-system integration.
- R console verbosity/warning scaffolding and expression-return helpers such as
  `r_pois_expr`.
- `.1` scalar versus vectorized duplicate entry points where the Fortran scalar
  routine already provides the computational kernel.

See `PORTING_NOTES.md` for historical-algorithm variants that intentionally
share a robust Fortran implementation.
