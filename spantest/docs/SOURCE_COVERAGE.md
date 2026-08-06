# Source coverage

## Exported R functions

- `R/span_bj.R` -> `src/spantest_classical.f90:span_bj`
- `R/span_f1.R` -> `src/spantest_classical.f90:span_f1`
- `R/span_f2.R` -> `src/spantest_classical.f90:span_f2`
- `R/span_grs.R` -> `src/spantest_classical.f90:span_grs`
- `R/span_hk.R` -> `src/spantest_classical.f90:span_hk`
- `R/span_km.R` -> `src/spantest_classical.f90:span_km`
- `R/span_py.R` -> `src/spantest_classical.f90:span_py`
- `R/span_gl_a.R` -> `src/spantest_gl.f90:span_gl_a`
- `R/span_gl_ad.R` -> `src/spantest_gl.f90:span_gl_ad`
- `R/span_as.R` -> `src/spantest_as.f90:span_as`
- `R/span_simulate.R` -> `src/spantest_simulation.f90:span_simulate`

## Internal computational helpers

- `f_ranklex` -> `rank_last_lex` in `spantest_gl` (private)
- `f_cauchypv` -> `cauchy_pvalue` in `spantest_as`
- `f_prods` -> `random_product_weights` in `spantest_as` (private)
- `f_ttest` / `f_testbm` -> `student_subseries_pvalues` and Cauchy merging in
  `spantest_as` (private)
- `f_getpv` / `f_getpv_batch` -> batched residual and score construction in
  `span_as`
- `f_rsstd` -> `standardized_skew_t` in `spantest_simulation`
- `src/gl_kernel.cpp:gl_sim_stats` -> streaming simulation loop in
  `spantest_gl`
- `src/sim_kernel.cpp:garch_filter` -> `garch_filter` in
  `spantest_simulation`

## R-only files omitted

- `R/RcppExports.R` and `src/RcppExports.cpp`: registration glue
- `R/spantest-package.R`: package metadata and documentation hooks
- `src/Makevars*`: R build-system configuration

The upstream package contains no plotting implementation.
