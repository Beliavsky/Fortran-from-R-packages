# Computational coverage

| Upstream routine | Fortran interface | Status |
|---|---|---|
| `edge` | `edge`, `edge_estimate` | Complete |
| `edge_rolling` | `edge_rolling` | Complete; fixed, adaptive, and endpoint widths |
| `edge_expanding` | `edge_expanding` | Complete |
| `spread` | `spread`, `spread_expanding`, `spread_endpoints` | Complete numerical dispatch |
| `EDGE` | `edge_estimate` and window routines | Complete |
| `AR`, `AR2` | `ar_estimate` | Complete |
| `CS`, `CS2` | `cs_estimate` | Complete |
| `ROLL` | `roll_estimate` | Complete |
| `OHLC` | `ohlc_estimate` | Complete, including dot-separated combinations |
| `sim` | `sim`, `simulate_ohlc` | Complete numerical simulation |
| `rmean`, `rsum`, `rfun` | Window slicing and NaN-aware statistics | Replaced by typed Fortran implementation |

## Excluded presentation infrastructure

- `data.table` rolling back end
- `xts` and `zoo` index preservation
- R data frames and column-name normalization
- R documentation and package registration at runtime

The original package remains available under `original/bidask-2.1.5`.
