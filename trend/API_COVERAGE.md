# API coverage

Status: **substantial**. **19 of 19 (100%)** exported computational R functions
have substantial numerical mappings. This is a function-mapping fraction, not
full R compatibility. Source inventory is the pinned local NAMESPACE, not a
cached website. All names are exact and match fpm.toml and README.md.

| R function | Public Fortran procedure | Mapping status |
|---|---|---|
| `bartels.test` | `trend_module:bartels_test` | substantial |
| `br.test` | `trend_module:br_test` | substantial |
| `bu.test` | `trend_module:bu_test` | substantial |
| `cs.test` | `trend_module:cs_test` | substantial |
| `csmk.test` | `trend_module:csmk_test` | substantial |
| `itam` | `trend_module:itam` | substantial |
| `lanzante.test` | `trend_module:lanzante_test` | substantial |
| `mk.test` | `trend_module:mk_test` | substantial |
| `mult.mk.test` | `trend_module:mult_mk_test` | substantial |
| `partial.cor.trend.test` | `trend_module:partial_cor_trend_test` | substantial |
| `partial.mk.test` | `trend_module:partial_mk_test` | substantial |
| `pettitt.test` | `trend_module:pettitt_test` | substantial |
| `rrod.test` | `trend_module:rrod_test` | substantial |
| `sea.sens.slope` | `trend_module:sea_sens_slope` | substantial |
| `sens.slope` | `trend_module:sens_slope` | substantial |
| `smk.test` | `trend_module:smk_test` | substantial |
| `snh.test` | `trend_module:snh_test` | substantial |
| `wm.test` | `trend_module:wm_test` | substantial |
| `ww.test` | `trend_module:ww_test` | substantial |

## Exclusions and counting

The 19 `export(...)` declarations are counted once each. `itam` performs numerical
slope/uncertainty calculations and is included; its registered plotting method
`plot.itam` is excluded. `rrod.test` dispatches to a **local** computational
method and is counted. `rrod.test.default` and `rrod.test.formula` are registered,
not independently exported, and are not additional entries. Plotting,
`summary.smktest`, datasets, and the internal `utilfn.R` helpers are outside the
coverage basis. No exported computational function is left untranslated.

See README.md for all interfaces, input minima, missing-value policies, statuses,
RNG differences, rank-bound behavior, stable-arithmetic corrections, and omitted
R metadata. Manifest mapping notes provide operation-specific distinctions.

## Reuse inspection

Before implementation the repository was checked for compatible shared APIs.
The implementation uses rfortran-core, Kendall, randtests, and exactRankTests
through sibling FPM paths; no dependency implementation is copied. Bartels beta
inference is compatible with randtests, while the normal approximation and
Cox-Stuart definition differ and retain trend formulas. Lanzante needs an exact
Wilcoxon branch absent from the existing stats translation. See NOTICE.md for
license compatibility and the extraDistr affine-beta substitution.
