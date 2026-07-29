# API reference

All public names are available through `use etrm`. Real calculations use
`real(dp)`, where `dp = kind(1.0d0)`.

Every main numerical routine returns an integer `status` and allocatable
`message`. A successful call has `status == etrm_ok`.

## Portfolio-insurance routines

### `cppi`

```fortran
call cppi(q, f, target_percent, risk_percent, result, status, message, &
          transaction_cost, integer_trades)
```

`risk_percent` may be a scalar or a vector matching `f`. As in the R source,
CPPI converts it to a currency risk factor using the initial futures price.

### `dppi`

```fortran
call dppi(q, f, target_percent, risk_percent, result, status, message, &
          transaction_cost, integer_trades)
```

A scalar risk percentage produces a constant currency risk factor based on the
initial price. A vector is multiplied elementwise by the futures-price vector.
The target moves only in the favorable direction.

### `obpi`

```fortran
call obpi(q, f, volatility, days_left, result, status, message, &
          strike_price, interest_rate, trading_days, &
          transaction_cost, integer_trades)
```

Defaults:

- `strike_price = f(1)`;
- `interest_rate = 0`;
- `trading_days = 250`;
- `transaction_cost = 0`; and
- `integer_trades = .true.`.

The routine uses Black-76 call deltas for buyers and put deltas for sellers.

### `shpi`

```fortran
call shpi(q, f, days_left, target_percent, result, status, message, &
          transaction_cost, integer_trades)
```

The initial schedule hedges `i/days_left` of the volume at observation `i`. Once
the portfolio reaches the cap or floor, the remaining position is fully hedged.

### `slpi`

```fortran
call slpi(q, f, target_percent, result, status, message, &
          transaction_cost, integer_trades)
```

The position remains unhedged until the cap or floor is reached, then becomes
fully hedged.

## `strategy_result`

Public fields:

- `name`, `volume`, `transaction_cost`, `integer_trades`;
- OBPI metadata: `strike_price`, `annual_volatility`, `interest_rate`,
  `trading_days`;
- arrays: `market`, `trade`, `exposed`, `position`, `hedge`, `target`,
  `portfolio`; and
- `risk_factor` for CPPI/DPPI.

The type-bound function `result%churn_rate()` returns
`sum(abs(trade))/abs(volume)`.

```fortran
call summarize_strategy(result, summary)
```

`strategy_summary` contains `volume`, `churn`, and seven-element vectors
`first`, `maximum`, `minimum`, and `last`. Their order is market, trade, exposed,
position, hedge, target, portfolio.

## Maximum-smoothness forward curve

Three generic forms are supplied:

```fortran
call msfc(include, contract, start_day, end_day, price, &
          result, status, message)

call msfc(include, contract, start_day, end_day, price, scalar_prior, &
          result, status, message)

call msfc(include, contract, start_day, end_day, price, prior_vector, &
          result, status, message)
```

`start_day` and `end_day` are integer offsets from the curve trade date. A vector
prior must cover day zero through the largest included `end_day`.

The descriptive core procedure is also public:

```fortran
call maximum_smoothness_forward_curve(include, contract, start_day, &
     end_day, price, prior_vector, result, status, message)
```

## `msfc_result`

- `n_days`, `n_contracts`, `n_polynomials`;
- `day`: offsets from zero through the final delivery end day;
- `curve`, `prior`;
- `knots` in years;
- `coefficients(5,n_polynomials)` for
  `a*t^4 + b*t^3 + c*t^2 + d*t + e`;
- selected `contract`, `original_index`, `start_day`, `end_day`;
- `market_price`; and
- `computed_price`, calculated by continuous integration of the curve over each
  delivery interval.

## Date helpers

```fortran
serial = civil_to_day(year, month, day)
offset = day_offset(year, month, day, base_year, base_month, base_day)
```

`civil_to_day` returns an `integer(int64)` day serial; `day_offset` returns a
default integer suitable for MSFC input.

## Status constants

- `etrm_ok`
- `etrm_err_size`
- `etrm_err_argument`
- `etrm_err_allocation`
- `etrm_err_linear_solve`
