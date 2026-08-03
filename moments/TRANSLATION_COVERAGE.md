# Translation coverage

Upstream namespace version: `moments` 0.14.1.

| R export | Fortran procedure | Status |
|---|---|---|
| `moment` | `moment` | Complete; vector and matrix overloads |
| `all.moments` | `all_moments` | Complete; vector and matrix overloads |
| `raw2central` | `raw2central` | Complete; vector and matrix overloads |
| `central2raw` | `central2raw` | Complete; vector and matrix overloads |
| `all.cumulants` | `all_cumulants` | Complete; corrected default plus legacy mode |
| `skewness` | `skewness` | Complete; vector and matrix overloads |
| `kurtosis` | `kurtosis` | Complete; vector and matrix overloads |
| `geary` | `geary` | Complete; vector and matrix overloads |
| `agostino.test` | `agostino_test` | Complete |
| `anscombe.test` | `anscombe_test` | Complete |
| `bonett.test` | `bonett_test` | Complete |
| `jarque.test` | `jarque_test` | Complete |

There is no plotting, compiled native code, formula interface, or external data
layer in the upstream package. All computational exports are included.
