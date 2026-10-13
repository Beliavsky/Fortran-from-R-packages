# API coverage

Package status: **substantial**. All **6 of 6 (100%)** distinct exported
computational R functions have numerical mappings. This fraction measures
mapped functions, not full R interface compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| alpha.proxy | burstfin_portfolio | alpha_proxy, alpha_proxy_grid | substantial |
| factor.model.stat | burstfin_factor | factor_model_stat, factor_covariance | substantial |
| slideWeight | burstfin_portfolio | slide_weight | complete |
| var.add.benchmark | burstfin_portfolio | var_add_benchmark | substantial |
| var.relative.benchmark | burstfin_portfolio | var_relative_benchmark | substantial |
| var.shrink.eqcor | burstfin_shrink | var_shrink_eqcor | substantial |

There are no unmapped functions in this basis. `partial.rainbow` is a color
presentation helper; `threeDarr` assembles arrays using R names/dimnames and
replication, an R-specific interface operation without a financial numerical
algorithm. Both are excluded, with their source retained for audit.
`fitted.statfacmodBurSt` is registered but not explicitly exported in NAMESPACE;
it is not counted separately. Its full/systematic/specific reconstruction is
implemented by `factor_covariance`.

The umbrella module `burstfin` re-exports the public operations and `dp`.
The complete argument contracts are documented next to each dummy declaration
with FORD comments; README.md describes typed outputs and departures from R.
