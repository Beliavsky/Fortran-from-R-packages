# Translation coverage

R names use dots; Fortran names use underscores and explicit argument arrays.
Closely related R front ends are intentionally consolidated when they differ
only by optional weights/household-size arguments.

## GB2 distribution, moments and indicators

| R export | Fortran API | Status |
| --- | --- | --- |
| `dgb2` | `dgb2` | translated |
| `pgb2` | `pgb2` | translated |
| `qgb2` | `qgb2` | translated |
| `rgb2` | `rgb2` | translated |
| `moment.gb2` | `moment_gb2` | translated |
| `incompl.gb2` | `incomplete_moment_gb2` | translated |
| `el.gb2` | `expected_log_gb2` | translated |
| `vl.gb2` | `variance_log_gb2` | translated |
| `sl.gb2` | `skewness_log_gb2` | translated |
| `kl.gb2` | `kurtosis_log_gb2` | translated |
| `arpt.gb2` | `arpt_gb2` | translated |
| `arpr.gb2` | `arpr_gb2` | translated |
| `rmpg.gb2` | `rmpg_gb2` | translated |
| `qsr.gb2` | `qsr_gb2` | translated |
| `Thomae` | `thomae_gb2` | translated |
| `gb2.gini`, `gini.gb2` | `gb2_gini`, `gini_gb2` | translated |
| `gini.b2` | `gini_b2` | translated |
| `gini.dag` | `gini_dag` | translated |
| `gini.sm` | `gini_sm` | translated |
| `main.gb2` | `main_gb2` | translated |
| `main2.gb2` | `main2_gb2` | translated |

## Likelihood and fitting

| R export | Fortran API | Status |
| --- | --- | --- |
| `logf.gb2` | `logf_gb2` | translated |
| `dlogf.gb2` | `dlogf_gb2` | translated |
| `d2logf.gb2` | `d2logf_gb2` | translated |
| `loglp.gb2`, `loglh.gb2` | `loglik_gb2(...,w,hs)` | consolidated/translated |
| `scoresp.gb2`, `scoresh.gb2` | `scores_gb2(...,w,hs)` | consolidated/translated |
| `info.gb2` | `info_gb2` | translated |
| `prof.gb2` | `prof_gb2` | translated |
| `proflogl.gb2` | `profile_loglik_gb2` | translated |
| `profscores.gb2` | `profile_scores_gb2` | translated |
| `fisk`, `fiskh` | `fisk_start(...,w,hs)` | consolidated/translated |
| `ml.gb2`, `mlh.gb2` | `fit_gb2_full(...,w,hs)` | consolidated/translated |
| `profml.gb2` | `fit_gb2_profile` | translated |
| `main.emp` | `main_emp` | translated without `laeken` |
| `mlfit.gb2` | `mlfit_gb2` | translated |
| `nlsfit.gb2` | `nlsfit_gb2` | translated |
| `robwts` | `robust_weights` | translated |

## Variance estimation

| R export | Fortran API | Status |
| --- | --- | --- |
| `varscore.gb2` | `varscore_gb2` | translated |
| `vepar.gb2` | `vepar_gb2` | translated |
| `derivind.gb2` | `derivind_gb2` | translated |
| `veind.gb2` | `veind_gb2` | translated |

Survey-design score variance is additionally exposed as
`survey_score_variance`, accepting a translated `survey_design_type`.

## Compound GB2

| R export | Fortran API | Status |
| --- | --- | --- |
| `fg.cgb2` | `fg_cgb2` | translated |
| `dl.cgb2` | `dl_cgb2` | translated |
| `pl.cgb2` | `pl_cgb2` | translated |
| `dcgb2` | `dcgb2` | translated |
| `pcgb2` | `pcgb2` | translated |
| `prcgb2` | `prcgb2` | translated |
| `mkl.cgb2` | `component_moments_cgb2` | translated |
| `moment.cgb2` | `moment_cgb2` | translated |
| `incompl.cgb2` | `incomplete_moment_cgb2` | translated |
| `qcgb2` | `qcgb2` | translated |
| `rcgb2` | `rcgb2` | translated |
| `arpt.cgb2` | `arpt_cgb2` | translated |
| `arpr.cgb2` | `arpr_cgb2` | translated |
| `rmpg.cgb2` | `rmpg_cgb2` | translated |
| `qsr.cgb2` | `qsr_cgb2` | translated |
| `main.cgb2` | `main_cgb2` | translated |
| `vofp.cgb2` | `vofp_cgb2` | translated |
| `pofv.cgb2` | `pofv_cgb2` | translated |
| `logl.cgb2` | `loglik_cgb2` | translated |
| `scoreU.cgb2` | `scoreu_cgb2` | translated |
| `scores.cgb2` | `scores_cgb2` | translated |
| `ml.cgb2` | `fit_cgb2` | translated |
| `hess.cgb2` | `hess_cgb2` | translated |
| `varscore.cgb2` | `varscore_mixture` | translated |
| `desvar.cgb2` | `survey_score_variance` | translated/consolidated |
| `vepar.cgb2` | `vepar_mixture` | translated |
| `derivind.cgb2` | `derivind_cgb2` | translated |
| `veind.cgb2` | `veind_cgb2` | translated |

## Compound GB2 with auxiliary information

| R export | Fortran API | Status |
| --- | --- | --- |
| `pkl.cavgb2` | `pkl_cavgb2` | translated |
| `lambda0.cavgb2` | `lambda0_cavgb2` | translated |
| `logl.cavgb2` | `loglik_cavgb2` | translated |
| `scores.cavgb2` | `scores_cavgb2` | translated |
| `scoreU.cavgb2` | `scoreu_cavgb2` | translated |
| `scorez.cavgb2` | `scorez_cavgb2` | translated |
| `ml.cavgb2` | `fit_cavgb2` | translated |
| `hess.cavgb2` | `hess_cavgb2` | translated |
| `varscore.cavgb2` | `varscore_mixture` / score helpers | translated/consolidated |
| `desvar.cavgb2` | `survey_score_variance` | translated/consolidated |
| `vepar.cavgb2` | `vepar_mixture` | translated/consolidated |
| `veind.cavgb2` | `veind_cavgb2_groups` | translated for group-dummy auxiliary design |

## Intentionally omitted plotting/presentation exports

- `plotsML.gb2`
- `saveplot`
- `contprof.gb2`
- `contindic.gb2`
- `dplot.cgb2`
- `dplot.cavgb2`

These routines create graphics or manage graphics devices and do not contribute
new statistical computations.
