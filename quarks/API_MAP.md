# API map

| R `quarks` function | Fortran counterpart | Coverage |
|---|---|---|
| `ewma` | `ewma` | Direct computational translation |
| `hs(method="plain")` | `hs(..., method_plain)` | Direct, including R type-7 quantiles and strict ES tail |
| `hs(method="age")` | `hs(..., method_age)` | Direct age weights, weighted interpolation, and tail weighting |
| `vwhs(model="EWMA")` | `vwhs(..., volatility_ewma)` | Direct |
| `vwhs(model="GARCH")` | `vwhs(..., volatility_garch)` | Gaussian sGARCH(1,1) through vendored `rugarch` |
| `fhs(model="EWMA")` | `fhs(..., volatility_ewma)` | Direct with explicit seedable RNG |
| `fhs(model="GARCH")` | `fhs(..., volatility_garch)` | Gaussian sGARCH(1,1) through vendored `rugarch` |
| `rollcast` | `rollcast` | Direct rolling alignment; smoothing adapted |
| `cvgtest` | `cvgtest` | Direct upstream formula by default; corrected transition option available |
| `trftest` | `trftest` | Direct binomial cumulative probability |
| `lossfun` | `lossfun` | Direct |
| `plop` | `plop`, `plop_time_varying` | Direct compatibility mode plus corrected exact-return mode |
| `plot.quarks`, `print.quarks` | none | R presentation infrastructure omitted |
| `runFTSdata` | none | Shiny/Yahoo download application omitted |

## Result types

- `risk_result`
- `rollcast_result`
- `coverage_result`
- `traffic_result`
- `loss_result`
- `pl_result`
- `rng_state`

Status constants such as `quarks_ok`, `quarks_empty_tail`, and
`quarks_no_violations` make edge behavior explicit without R exceptions.
