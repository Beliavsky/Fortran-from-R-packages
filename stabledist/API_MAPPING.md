# API mapping

| Upstream R | Fortran | Notes |
|---|---|---|
| `dstable` | `dstable` | Scalar and vector overloads; S0/S1/S2; log density |
| `pstable` | `pstable` | Scalar/vector; lower/upper tail; log probability |
| `qstable` | `qstable` | Scalar/vector; lower/upper tail; log probability input |
| `rstable` | `rstable` | Scalar scale/location parameters |
| `rstable` with vector `gamma`,`delta` | `rstable_varying` | Explicit recycled-vector API |
| `stableMode` | `stable_mode` | Bounded numerical maximization of S0 density |
| `.om` | `stable_omega` | S0/S1 parameter conversion helper |
| `C.stable.tail` | `c_stable_tail` | Scalar/vector overloads; optional log form |
| `dPareto` | `dpareto` | Upstream tail-density approximation |
| `pPareto` | `ppareto` | Positive-tail helper; upstream negative-x FIXME is not emulated |
| `.fct1`, `.fct2` | internal | Nolan density integral kernels |
| `.FCT1`, `.FCT2` | internal | Nolan CDF integral kernels |
| `.integrate2` | internal `integrate_split` + `r_mod::integrate` | No R warning object; numerical value only |
| `.unirootNA` | internal `bisect_root` | Bracketed scalar root solver |

R options/debug printing, S3 behavior, plotting examples, and package-check
helpers such as `doExtras` are intentionally omitted.
