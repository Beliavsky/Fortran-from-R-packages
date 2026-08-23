# Translation coverage

| Upstream R export | Fortran API | Status |
|---|---|---|
| Runif | `runif` | translated |
| Rbeta | `rbeta` | translated |
| Rexp | `rexp` | translated |
| Rchisq | `rchisq` | translated |
| Rgamma | `rgamma` | translated |
| Rgeom | `rgeom` | translated |
| Rcauchy | `rcauchy` | translated; upstream bug corrected |
| Rt | `rt` | translated |
| Rpareto | `rpareto` | translated |
| Rfrechet | `rfrechet` | translated |
| Rlaplace | `rlaplace` | translated |
| Rgumble | `rgumble`, `rgumbel` | translated |
| Rarcsine | `rarcsine` | translated |
| Rnorm | `rnorm` | translated; Box-Muller replaces external zigg backend |
| all `*.mat` generators | corresponding `*_mat` | translated |
| all `colR*` generators | corresponding `col_r*` | translated |
| Sample.int | `sample_int` | translated; bounds bugs corrected |
| Sample | `sample_real` | translated |
| colSample | `col_sample` | translated numerically; serial execution |
| rowSample | `row_sample` | translated numerically; serial execution |
| nanoTime | `nano_time` | translated as portable elapsed-clock ticks |
| setSeed | `set_seed` | translated with deterministic PCG streams |

There is no plotting or statistical-model layer in the upstream package.
