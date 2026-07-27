# Computational coverage

All exported numerical functions in `FinancialMath` 0.1.1 are represented.
Names containing dots in R use underscores in Fortran.

| Upstream R function | Fortran API | Notes |
|---|---|---|
| `NPV` | `npv` | Arbitrary positive times |
| `IRR` | `irr` | Bracketed real root search |
| `TVM` | `solve_tvm` | Explicit unknown selector |
| `cf.analysis` | `cf_analysis` | PV, Macaulay/modified duration and convexity |
| `rate.conv` | `rate_conv` | Interest, discount, and force conversion |
| `yield.dollar` | `yield_dollar` | Dollar-weighted yield |
| `yield.time` | `yield_time` | Time-weighted linked yield |
| `swap.rate` | `swap_rate` | Spot-rate or ZCB-price input |
| `swap.commodity` | `swap_commodity` | Discount-weighted fixed commodity price |
| `annuity.level` | `annuity_level` | PV/FV/payment/period/rate solving |
| `annuity.arith` | `annuity_arith` | Arithmetic-gradient annuity |
| `annuity.geo` | `annuity_geo` | Geometrically growing annuity |
| `perpetuity.level` | `perpetuity_level` | Immediate or due |
| `perpetuity.arith` | `perpetuity_arith` | Arithmetic gradient |
| `perpetuity.geo` | `perpetuity_geo` | Growing perpetuity |
| `amort.table` | `amort_table` | Typed complete schedule |
| `amort.period` | `amort_period` | Payment decomposition at a period |
| `bond` | `bond` | Price, premium/discount, durations, convexities |
| `bls.order1` | `bls_order1` | Call/put price, delta, theta, and vega |
| `option.call` | `option_call` | Long/short payoff and profit table |
| `option.put` | `option_put` | Long/short payoff and profit table |
| `covered.call` | `covered_call` | Numerical payoff and profit arrays |
| `covered.put` | `covered_put` | Upstream payoff convention retained |
| `protective.put` | `protective_put` | Numerical payoff and profit arrays |
| `straddle` | `straddle` | Caller-supplied premiums |
| `straddle.bls` | `straddle_bls` | Black-Scholes premiums |
| `strangle` | `strangle` | Caller-supplied premiums |
| `strangle.bls` | `strangle_bls` | Black-Scholes premiums |
| `bull.call` | `bull_call` | Caller-supplied premiums |
| `bull.call.bls` | `bull_call_bls` | Black-Scholes premiums |
| `bear.call` | `bear_call` | Caller-supplied premiums |
| `bear.call.bls` | `bear_call_bls` | Black-Scholes premiums |
| `collar` | `collar` | Caller-supplied premiums |
| `collar.bls` | `collar_bls` | Black-Scholes premiums |
| `butterfly.spread` | `butterfly_spread` | Caller-supplied premiums |
| `butterfly.spread.bls` | `butterfly_spread_bls` | Black-Scholes premiums |
| `forward` | `forward_contract`, `forward_price` | None, continuous, or discrete dividends |
| `forward.prepaid` | `forward_prepaid`, `prepaid_forward_price` | Discounted delivery price |

## Excluded presentation-only behavior

The R package's plotting branches, device management, labels, legends, and
R data-frame/matrix row names are not compiled. Every underlying value used by
those plots is available through a scalar result or typed array.
