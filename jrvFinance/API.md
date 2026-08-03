# API

Import the aggregate module:

```fortran
use jrvfinance
```

## Core types

- `date_t`: Gregorian year, month, and day
- `root_result`: root, residual, iteration count, status, and message
- `bond_cashflows`: payment times, cash flows, accrued interest, and status
- `annuity_breakup_result`: opening/closing principal and payment components
- `black_scholes_result`: call/put values, eight Greeks, and intermediate values

Status constants include `JRV_OK`, `JRV_INVALID_ARGUMENT`,
`JRV_NO_CONVERGENCE`, `JRV_NO_ROOT`, `JRV_NONFINITE`, and
`JRV_SIZE_MISMATCH`.

## Discounted cash flow

- `equiv_rate(rate, from_freq, to_freq)`
- `npv(cf, rate, cf_freq, comp_freq, cf_t, immediate_start, status)`
- `irr(cf, interval, cf_freq, comp_freq, cf_t, r_guess, toler, convergence, max_iter, method, status)`
- `duration(cf, rate, cf_freq, comp_freq, cf_t, immediate_start, modified, status)`

Use frequency `0.0_dp` for continuous compounding.

## Annuities

- `annuity_pv`
- `annuity_fv`
- `annuity_instalment`
- `annuity_periods`
- `annuity_rate`
- `annuity_instalment_breakup`

Periods are real-valued to support fractional terms and an IEEE positive
infinity default for perpetuities.

## Dates and bonds

- `date`, `parse_date`, `date_string`, `edate`
- `daycount_actual`, `daycount_30_360`, `year_fraction`
- `coupons_n`, `coupons_next`, `coupons_prev`, `coupons_dates`
- `bond_tcf`, `bond_price`, `bond_yield`, `bond_duration`
- `bond_prices`, `bond_yields`, `bond_durations`

Supported bond conventions are `30/360`, `30/360E`, `ACT/360`, and `ACT/ACT`.
Coupon frequencies must divide 12.

## Options

- `gen_bs(s, strike, rate, sigma, time, div_yield)`
- `gen_bs_implied(s, strike, rate, price, time, div_yield, put_option, ...)`
- `gen_bs_implied_guess` (exposed helper)

`gen_bs` returns call and put prices, deltas, thetas, gamma, vega, rhos,
exercise probabilities, and the usual `d1`/`d2` quantities.

## Root solvers

Callbacks implement:

```fortran
subroutine callback(x, context, value, gradient)
  real(dp), intent(in) :: x
  class(*), intent(in) :: context
  real(dp), intent(out) :: value, gradient
end subroutine
```

The public solvers are `newton_raphson_root`, `bisection_root`, and
`irr_solve`.
