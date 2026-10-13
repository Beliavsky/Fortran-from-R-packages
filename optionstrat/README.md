# optionstrat: European option and spread analytics

Modern free-form Fortran translation of CRAN **optionstrat 1.4.1**, by John T. Buynak.
GPL-3.0-only. See [NOTICE.md](NOTICE.md), [LICENSE](LICENSE), and the preserved
[upstream source](upstream/R/optionstrat.R). This is a numerical library, not investment advice.

## Build and use

Extract as the top-level directory `optionstrat` in
[Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages), alongside
`rfortran-core`. With gfortran and FPM on Windows, Linux, or macOS:

```sh
fpm build
fpm test
fpm run --example option_analysis
```

The only dependency is `rfortran-core = { path = "../rfortran-core" }`.
Its `r_kinds:dp`, normal density/CDF/quantile functions, and stable `r_log1p` are reused;
no dependency source, BLAS/LAPACK library, R installation, or system link flags are needed.
The exported `dp` is the same kind as the dependency. Do not enable `-ffast-math`,
`-Ofast`, `-ffinite-math-only`, or equivalent flags: NaN/infinity behavior is part of the API.

```fortran
use optionstrat, only : dp, option_result, calleval, implied_volatility
! Declarations belong before executable statements.
type(option_result) :: opt
real(dp) :: sigma
opt = calleval(100.0_dp, 100.0_dp, 0.2_dp, 0.25_dp, 0.03_dp, 0.01_dp)
sigma = implied_volatility('call', opt%premium, 100.0_dp, 100.0_dp, 0.25_dp, 0.03_dp, 0.01_dp)
```

## Important source conventions and safer alternatives

This translation deliberately distinguishes compatibility with the source from corrected calculations:

- `calltheta` and `puttheta`, including `theta` fields and Greek selectors, preserve the
  upstream formula. Its interest and dividend terms are divided by tenor in addition to
  the daily conversion, so it is generally **not** the price's standard daily time derivative.
  Use `calltheta_standard` and `puttheta_standard` for conventional daily theta.
- `iv_calc` preserves the source nearest-price search on volatility `0, 0.001, ..., 5`.
  It does not enforce arbitrage bounds: an out-of-range premium selects a boundary grid value.
  If forward log-moneyness is exactly zero, the R zero-volatility candidate is NaN and R's
  `min`/`match` path chooses the first grid point. `iv_calc` preserves the resulting **zero IV**,
  even for a positive premium. Use `implied_volatility` for accurate bounded inversion.
- `prob_above_quantile` negates the **whole** normal quantile, including its mean.
  `prob_btwn_quantile` uses the negative of its lower change as its upper change. This is
  the source convention and gives the wrong mean shift when `mean` is nonzero.
  Use `prob_above_quantile_standard` and `prob_btwn_quantile_standard` for corrected quantiles.
- Spread probabilities use the source's normal distribution of arithmetic returns, not the
  lognormal risk-neutral distribution of the option-pricing model. These are distinct models.

The corrected extensions are not extra R-function mappings. Tests compare compatibility
operations against genuine R and independently verify the corrected derivatives/inversions.

## Units, formulas, and inputs

Option routines use positional `(s, x, sigma, t, r, d)`, with `d` optional and zero by default.
Spot and strike are positive currency prices. Volatility, interest, and dividend inputs are
**annual decimals**: 20% volatility is `0.20_dp`; 2% continuously compounded rate is `0.02_dp`.
Time is years. Rho and vega are price changes per **one percentage point**, so each is
`0.01` times the derivative with respect to the corresponding decimal input. Delta is price
change per unit spot; gamma is delta change per unit spot. Standard theta is per calendar day,
using 365 days per year. R-compatible theta has the extra tenor scaling described above.

For positive volatility, `d1 = [log(s/x) + (r-d)t] / (sigma sqrt(t)) + sigma sqrt(t)/2`,
`d2 = d1 - sigma sqrt(t)`, and the price is the discounted standard Black-Scholes-Merton formula.
Logarithms are differenced internally to avoid overflow in `s/x`. The source's `t == 0`
replacement by `1e-200` is preserved, including extremely large ATM gamma/theta. Zero time
is therefore a tiny-time convention, not a definition of all derivatives at expiry.

Zero volatility is an explicit extension: discounted intrinsic value; delta equals the
appropriate discounted indicator, or its midpoint at the deterministic forward-price kink;
gamma is zero away from the kink and NaN at it. Vega uses the positive-volatility right limit.
The same rule applies when `sigma * sqrt(t)` underflows to zero. Undefined zero-price
elasticity returns NaN instead of R's possible infinities. Negative rates/dividend yields
are accepted if discounts are representable; invalid signs, nonfinite inputs, and detected
discount-factor overflow return NaN. At extreme representability limits some derivatives
may overflow. No warnings, exceptions, or R coercions are emitted.

`implied_volatility(type, price, s, x, t, r, d)` uses an expanding bracket then bisection.
It requires positive tenor and a premium in the no-arbitrage interval [discounted intrinsic,
discounted spot/strike). The lower bound returns zero; the upper bound, invalid type, and
out-of-bounds premiums return NaN. At a price rounded exactly to the lower bound volatility
is numerically unidentifiable, and zero is returned. `iv_calc` deliberately has different
compatibility behavior. Option types and Greek selectors are exact lowercase strings.

## API shapes and results

All public procedures are pure elemental, so scalar arguments can accompany conformable
arrays of any rank; derived-type results have the same shape. Empty arrays are supported.
Array arguments must conform under normal Fortran rules; R recycling is not implemented.

- `callpremium`, `putpremium`, `calldelta`, `putdelta`, `optiongamma`, `optionvega`,
  `calltheta`, `puttheta`, `callrho`, `putrho` return scalars or elemental arrays.
- `calleval` / `puteval` return `option_result`: `premium, delta, gamma, vega, theta, rho`.
  `opteval` returns `dual_option_result` with `call` and `put` components.
- `callgreek(greek, ...)` / `putgreek(greek, ...)` select `premium`, `delta`, `gamma`,
  `vega`, `theta`, or `rho`; unknown names return NaN. `lambda(type, ...)` gives signed elasticity.
- `prob_below(spot, upper, asd, dsd, dte, mean, tradedays)`,
  `prob_above(spot, lower, ...)`, and `prob_btwn(spot, lower, upper, ...)` return probabilities.
  Arguments after the threshold(s) are optional: zero except `tradedays = 262`.
  All volatilities and arithmetic-return changes are decimals, even though R labels them
  "percent change". `mean` is used as the horizon return location without rescaling by days.
  If only `asd > 0`, horizon SD is `asd / sqrt(tradedays) * sqrt(dte)`; if only `dsd > 0`
  and `dte > 0`, it is `dsd * sqrt(dte)`. If both are positive, the source uses `asd`
  directly as horizon SD. If `dte == 0` with daily-only volatility, the source uses zero SD.
  All-zero volatility/day inputs return NaN; zero SD otherwise uses a point-mass CDF,
  including the boundary in the lower tail. Reversed between bounds give a signed difference.
- Quantile counterparts take `(spot, p, asd, dsd, dte, mean, tradedays)`.
  Above/below routines return the price; their decimal change is `price / spot - 1`.
  Between routines return `price_range`: `probability, lower_change, upper_change,
  lower_price, upper_price`. Endpoints `p = 0, 1` can return infinities for positive SD.
- `vertical(options, s, x1, x2, t, r, sigma, sigma2, vol, d)` uses short strike `x1`
  and long strike `x2`. Optional `sigma2` and probability `vol` default to `sigma`, `d` to zero.
  `vertical_result` fields: `spot, short_strike, long_strike, max_profit, max_loss,
  breakeven, prob_be, prob_max_profit, prob_max_loss, initial_dc`.
  Maximum loss is signed P/L, usually negative; `initial_dc` is positive for a credit.
  The four source branches and breakeven conventions are retained without rounding.
  Equal strikes are rejected. Unequal leg volatilities can produce unusual economic results;
  the library preserves the formulas rather than imposing cross-strike arbitrage consistency.
- `dv(s, x1, x2, x3, x4, t, r, sigma, sigma2, sigma3, sigma4, vol, d)` represents a double
  vertical credit spread, requiring `x1 < x2 <= x3 < x4`. Extra volatilities default to `sigma`.
  `double_vertical_result` contains `spot, long_put, short_put, short_call, long_call,
  lower_be, higher_be, prob_be, prob_max_profit, prob_max_loss`. The original does not return
  maximum profit/loss amounts, and neither does this result. Both spread routines convert
  annual volatility with **256** trading days (`vol / 16`) and tenor with **365** calendar days.
- `r_cont(r, n)` uses the stable equivalent `n * log1p(r/n)`; `n > 0`, `1 + r/n > 0`.
  `tdiff(date1, date2, period)` accepts ISO Gregorian `YYYY-MM-DD` dates, years 0001–9999,
  and exact `days` or `years`. Years are day differences divided by 365, including leap days.
  Time zones, time of day, locale parsing, and arbitrary R date classes are omitted.

## Verification

`fpm test` runs 368 deterministic checks. Checked-in numerical fixtures were produced by
actual R optionstrat 1.4.1 using [test/generate_r_fixtures.R](test/generate_r_fixtures.R).
The embedded [reference module](test/r_reference_data.f90) needs no R at test time.
Fixtures cover every mapped operation, both option types and all Greeks, six regular and
two expiry scenarios, every vertical-spread branch, double verticals, probability scaling,
nonzero means, and quantile conventions. R spreads are rounded to two decimals, so comparison
allows half a cent; other comparisons use floating-point tolerances. Dates/rate conversions
are checked against exact independent values rather than an R runtime dependency.

Independent tests check a direct ERFC-derived Black-Scholes price, signed Greeks, put-call
parity, finite-difference delta/gamma/vega/rho and corrected theta, continuous-IV round trips,
quantile inversion, zero tenor/volatility, deterministic kinks, extreme moneyness, NaN/Inf and
bad-input handling, arbitrage bounds, Gregorian dates, tiny rates, and vector/matrix/empty shapes.
The source quirks are regression tested explicitly. Default/debug/release builds, tests, and
the example are part of release validation. No stochastic tests are used.

## Translation coverage

Package status: **substantial**. **24 of 24 (100%)** exported computational R functions are mapped.
This fraction measures mapped computational functions, not complete R compatibility.
The six plotting exports are excluded from the denominator. No computational export is left
untranslated. R graphics, data frames, transposition/rounding presentation, warnings/errors,
recycling/coercion, general date interfaces, and interactive/plot interfaces are omitted.
See [API_COVERAGE.md](API_COVERAGE.md) for function-specific differences.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `calldelta` | `optionstrat` | `calldelta` | substantial |
| `calleval` | `optionstrat` | `calleval` | substantial |
| `callgreek` | `optionstrat` | `callgreek` | substantial |
| `callpremium` | `optionstrat` | `callpremium` | substantial |
| `callrho` | `optionstrat` | `callrho` | substantial |
| `calltheta` | `optionstrat` | `calltheta` | substantial |
| `dv` | `optionstrat` | `dv` | substantial |
| `iv.calc` | `optionstrat` | `iv_calc` | substantial |
| `lambda` | `optionstrat` | `lambda` | substantial |
| `opteval` | `optionstrat` | `opteval` | substantial |
| `optiongamma` | `optionstrat` | `optiongamma` | substantial |
| `optionvega` | `optionstrat` | `optionvega` | substantial |
| `prob.above` | `optionstrat` | `prob_above`, `prob_above_quantile` | substantial |
| `prob.below` | `optionstrat` | `prob_below`, `prob_below_quantile` | substantial |
| `prob.btwn` | `optionstrat` | `prob_btwn`, `prob_btwn_quantile` | substantial |
| `putdelta` | `optionstrat` | `putdelta` | substantial |
| `puteval` | `optionstrat` | `puteval` | substantial |
| `putgreek` | `optionstrat` | `putgreek` | substantial |
| `putpremium` | `optionstrat` | `putpremium` | substantial |
| `putrho` | `optionstrat` | `putrho` | substantial |
| `puttheta` | `optionstrat` | `puttheta` | substantial |
| `r.cont` | `optionstrat` | `r_cont` | substantial |
| `tdiff` | `optionstrat` | `tdiff` | substantial |
| `vertical` | `optionstrat` | `vertical` | substantial |
