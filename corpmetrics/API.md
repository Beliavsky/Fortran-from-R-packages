# API

All public procedures are in module `corpmetrics`. Real calculations use
`dp = kind(1.0d0)`.

## High-level procedures

### `balsh`

```fortran
call balsh(fa, ca, inv, fl, cl, result, status)
```

Returns total assets, liabilities, equity, working capital, current ratio,
acid-test ratio, leverage, and debt-to-equity ratio in
`balance_sheet_result`.

### `capm`

```fortran
call capm(rf, ri, rm, result, status)
```

Computes sample beta, expected returns, CAPM required return, and valuation
classification in `capm_result`. `ri` and `rm` must have equal lengths of at
least two. `rf` may have any nonzero length, matching the executable R code.

Valuation constants:

- `capm_undervalued`
- `capm_fairly_valued`
- `capm_overvalued`

`valuation_label` converts a code to the original R label.

### `ddm`

```fortran
call ddm(dividend, required_return, result, g1, g2, period, status)
```

Valid optional-argument combinations are:

- no optionals: zero-growth model
- `g1`: Gordon model
- `g1`, `g2`, and `period`: differential-growth model

The result model constants are `ddm_zero_growth`, `ddm_gordon`, and
`ddm_differential`. `ddm_model_label` returns the original display label.

### `fis`

```fortran
call fis(face_value, coupon_rate, ytm, maturity, result, semiannual, status)
```

Returns price, Macaulay duration, modified duration, and payment-period count in
`fixed_income_result`. `maturity` is an integer number of years.

### `idm`

```fortran
call idm(cash_flows, costs, result, status)
```

Returns NPV and constant-rate IRR in `investment_result`. The NPV uses the
original element-specific costs:

```text
sum(cash_flows(t) / (1 + costs(t))^(t-1))
```

The IRR solves the same cash-flow stream using one common rate on the original
interval `(-1, 1]`.

Lower-level procedures:

```fortran
value = net_present_value(cash_flows, costs)
call internal_rate_of_return(cash_flows, irr, status, lower, upper, tolerance, max_iterations)
```

### `insta`

```fortran
call insta(revenue, cost_of_sales, net_income, result, &
    preferred_dividend, shares, price_per_share, status)
```

Returns gross and net margins. EPS is calculated only when both
`preferred_dividend` and `shares` are present. P/E is calculated only when EPS
and `price_per_share` are available. Flags `has_eps` and `has_pe` distinguish
missing optional results from numerical zero.

### `loan`

```fortran
call loan(amount, rate, periods, result, status)
```

Returns the unrounded installment and total repayment plus allocated arrays for
period, interest, principal, and balance. Table values are rounded to two
places after each period, matching the R presentation behavior; the recursive
balance itself uses unrounded arithmetic.

## Status codes

- `cm_success`
- `cm_invalid_input`
- `cm_dimension_mismatch`
- `cm_singular`
- `cm_root_not_bracketed`
- `cm_root_failure`

When a ratio or root is undefined, the result is set to IEEE NaN where useful.
