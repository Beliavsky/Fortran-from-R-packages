# API mapping

| Upstream relsurv routine/source | Fortran implementation |
|---|---|
| `pystep.c` | `pystep` |
| `pystep2.c` | `pystep2` |
| `exps.c`, `popsurv` | `expected_survival`, `population_survival_curve` |
| `expprep2` | `expprep2_expected`, `expprep2_summary` |
| `netwei.c` | `netwei_summary` |
| `netfastpinter2.c` | `netfast_summary` |
| `rs.surv` | `rs_surv` |
| `cmp.rel` / `cmpfast.c` | `cmp_rel` |
| `rs.diff` | `rsdiff` |
| `nessie` | `nessie_expected` |
| `aalen_beta.cpp` | `prepare_x`, `fit_ols2`, Aalen fit routines |
| `rsaalen` | `aalen_fit_relative` |
| `survaalen` | `aalen_fit` |
| `rsadd` direct ML | `rsadd_ml_rows`, `rsadd_piecewise` |
| `rsadd` EM | `rsadd_em`, `rsadd_em_core` |
| `rsadd` grouped binomial GLM | `rsadd_glm_bin` |
| `rsadd` grouped Poisson GLM | `rsadd_glm_poisson` |
| `residuals.rsadd` numerical core | `rsadd_schoenfeld_residuals` |
| `rs.br` | `rs_br` |
| `rs.zph` | `rs_zph` |
| `rstrans` | `rstrans_times`, `rstrans_fit` |
| `rsmul` | `rsmul_fit` |
| `survsplit` | `survsplit_counting` |
| `invtime` | `inverse_time_monotone` |
| `transrate` | `transrate` |
| `transrate.hld` | `transrate_hld` |
| `transrate.hmd` | `transrate_hmd` |
| `joinrate` | `join_ratetables` |
| `epa` | `epa_smooth`, `epanechnikov_boundary_matrix` |
| `years` default YD | `years_difference` |
| `years` YL2013 | `years_yl2013` |
| `years` YL2017 | `years_yl2017` |
| `years` Greenwood area variance | `greenwood_area_variance` |
| `years` bootstrap aggregation | `bootstrap_column_variance` plus optional replicate arrays |
