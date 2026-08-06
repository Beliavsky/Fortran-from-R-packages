# API map

| Upstream R routine | Fortran routine | Status |
|---|---|---|
| `mod_duration(yields, maturity)` | elemental `mod_duration(yield_rate, maturity, source_compatible)` | Implemented |
| `convexity(yields, maturity)` | elemental `convexity(yield_rate, maturity, source_compatible)` | Implemented |
| `total_return(yields, maturity, mdur, convex, scale)` | generic `total_return` for vectors and matrices | Implemented |
| internal one-period expression | elemental `period_total_return` | Added direct API |
| `get_yields(..., na_locf, percent_adjust)` | `prepare_yields`, `carry_forward`, `percent_to_decimal` | Computational preprocessing implemented; FRED download omitted |
| `xts_to_tibble` | none | Omitted: R container conversion |
| `tibble_to_xts` | none | Omitted: R container conversion |

## Data representation

R `xts` vectors are represented by rank-1 `real(dp)` arrays. Multiple `xts`
columns are represented by rank-2 arrays with shape `(time, series)`. Dates and
column names remain the responsibility of the calling application.

## Error handling

`total_return` optionally returns a status code and message:

- `tt_success`
- `tt_err_size`
- `tt_err_maturity`
- `tt_err_scale`

Invalid individual yield values produce IEEE quiet NaNs, matching the natural
propagation behavior of the R arithmetic.
