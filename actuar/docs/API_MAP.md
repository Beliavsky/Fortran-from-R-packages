# API map

Representative upstream-to-Fortran mappings in v0.3.0:

| R / C API | Fortran API |
|---|---|
| heavy-tail `d/p/q/r/m/lev` functions | corresponding scalar/elemental Fortran procedures |
| `dpoisinvgauss` and related | same names; upstream recurrence retained |
| `dphtype`, `pphtype`, `rphtype`, moments/MGF | same names |
| `panjer` | `panjer_ab` plus frequency-specific wrappers |
| `exact` aggregate method | `aggregate_exact` |
| normal / normal-power aggregate approximation | `aggregate_normal_*`, `aggregate_npower_*` |
| `rcompound` | `compound_sums`, `rcompound_callbacks` |
| `ruin` | `ruin_phase_type` |
| `bstraub` | `bstraub_fit` |
| `hierarc` basic matrix-first recursion | `hierarchical_credibility` |
| `hierarc` + `actuar_do_hierarc` exact recursion | `hierarc_exact_fit` |
| Hachemeister origin method | `hachemeister_fit` |
| Hachemeister barycenter method | `hachemeister_barycenter_fit` |
| `mde(..., measure="CvM")` | `mde_cvm`, `mde_grouped_cvm` |
| `mde(..., measure="chi-square")` | `mde_grouped_chisq` |
| `mde(..., measure="LAS")` | `mde_grouped_las` |
| `coverage()` CDF helper | `coverage_cdf` + `coverage_spec_t` |
| `coverage()` PDF/mass helper | `coverage_pdf` + `coverage_spec_t` |
| `rcomphierarc` / `simul` | `rcomphierarc_simulate` |
| `emm` | `empirical_moments`, `grouped_moments` |
| `ogive` | `ogive_eval` |
| `elev` | `elev_individual`, `elev_grouped` |
| `VaR.aggregateDist` | `aggregate_var` / type-bound `quantile` |
| `CTE.aggregateDist` | `aggregate_cte` |

R vector recycling, S3 classes, formulas and expression evaluation are not
replicated. Fortran arrays, derived result types and explicit callbacks are
used instead.
