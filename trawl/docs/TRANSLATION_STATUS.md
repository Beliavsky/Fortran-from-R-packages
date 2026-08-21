# Translation status

Upstream: `trawl` 0.2.2 (GPL-3), CRAN publication 2021-02-22.

## Export map

| R export | Fortran API | Status |
|---|---|---|
| `trawl_Exp` | `trawl_exp` | translated |
| `trawl_DExp` | `trawl_dexp` | translated |
| `trawl_supIG` | `trawl_supig` | translated |
| `trawl_LM` | `trawl_lm` | translated |
| `acf_Exp` | `acf_exp` | translated |
| `acf_DExp` | `acf_dexp` | translated |
| `acf_supIG` | `acf_supig` | translated |
| `acf_LM` | `acf_lm` | translated |
| `fit_Exptrawl` | `fit_exptrawl` | translated |
| `fit_supIGtrawl` | `fit_supigtrawl` | translated |
| `fit_LMtrawl` | `fit_lmtrawl` | translated |
| `fit_DExptrawl` | `fit_dexptrawl` | translated |
| `fit_marginalPoisson` | `fit_marginal_poisson` | translated |
| `fit_marginalNB` | `fit_marginal_nb` | translated |
| `fit_trawl_intersection` | `fit_trawl_intersection` | translated |
| `fit_trawl_intersection_Exp` | `fit_trawl_intersection_exp` | translated |
| `fit_trawl_intersection_LM` | `fit_trawl_intersection_lm` | translated |
| `Bivariate_NBsim` | `bivariate_nbsim` | translated |
| `Bivariate_LSDsim` | `bivariate_lsdsim` | translated |
| `Trivariate_LSDsim` | `trivariate_lsdsim` | translated |
| `LSD_Mean` | `lsd_mean` | translated |
| `LSD_Var` | `lsd_var` | translated |
| `ModLSD_Mean` | `modlsd_mean` | translated |
| `ModLSD_Var` | `modlsd_var` | translated |
| `BivLSD_Cor` | `bivlsd_cor` | translated |
| `BivLSD_Cov` | `bivlsd_cov` | translated |
| `BivModLSD_Cov` | `bivmodlsd_cov` | translated |
| `BivModLSD_Cor` | `bivmodlsd_cor` | translated |
| `sim_UnivariateTrawl` | `sim_univariate_trawl` | translated |
| `sim_BivariateTrawl` | `sim_bivariate_trawl` | translated |
| `plot_2and1hist` | - | omitted: plotting |
| `plot_2and1hist_gg` | - | omitted: plotting |

## Parity notes

The package is computationally self-contained rather than linking the R-package
imports. Statistical formulas and parameterizations are translated directly.
Random streams therefore do not reproduce R's seed-for-seed output.

The upstream `fit_DExptrawl` assigns `Delta <- 1` inside its objective and thus
ignores a non-default `Delta` during estimation. The Fortran routine preserves
this behavior by default for parity. Pass `preserve_upstream_delta_bug=.false.`
to use the supplied `delta_t` in the GMM objective.

The three GMM fitting routines now use the supplied `DEoptim-fortran` 0.1.0
translation of DEoptim 2.2-8 rather than the earlier simplified local
DE/rand/1/bin substitute. The wrapper reproduces upstream trawl's DEoptim
control values (`strategy=2`, `NP=10*npar`, `CR=0.5`, `F=0.8`, and 1000
iterations by default). Random streams are still not seed-for-seed identical to
R because the translated DEoptim engine uses a standalone RNG.

Trawl intersections use the same root search interval `[-1000,-1e-6]` and the
same 0/1/2/>=3-root decision logic. Segment integrals are evaluated from closed
form antiderivatives of the four trawl families instead of calling R's numerical
`integrate`, improving reproducibility without changing the mathematical target.
