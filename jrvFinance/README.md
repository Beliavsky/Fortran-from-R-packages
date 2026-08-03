# jrvFinance-fortran

A modern Fortran 2018 translation of the computational core of the R package
`jrvFinance` 1.4.3. The project is self-contained and builds with the Fortran
Package Manager (FPM).

## Implemented functionality

- NPV, IRR, equivalent-rate conversion, duration, and modified duration
- Present/future value, instalment, term, rate, and amortization breakup for annuities
- Gregorian dates, month shifts, actual and 30/360 day counts, and year fractions
- Coupon schedules, accrued interest, bond cash flows, prices, yields, and durations
- Generalized Black-Scholes call/put prices, Greeks, and implied volatility
- Safeguarded Newton-Raphson, geometric bracketing, and bisection root solvers

All 30 functions exported by the upstream R namespace have computational
Fortran equivalents. Dots in R names are replaced by underscores, for example
`bond.price` becomes `bond_price` and `annuity.pv` becomes `annuity_pv`.

## Build

```text
fpm build
fpm test
fpm run
```

The demo is `demo_jrvfinance`; examples are also available as FPM targets.

## API conventions

Dates use the `date_t` derived type and can be created with `date(year,month,day)`
or `parse_date("YYYY-MM-DD")`. Continuous compounding is represented by a
frequency of `0.0_dp`; positive values represent compounding periods per year.
R lists are represented by derived result types such as `black_scholes_result`,
`bond_cashflows`, and `root_result`.

See `API.md`, `PORTING.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

GPL-2.0-or-later, matching the upstream package declaration `GPL (>= 2)`.
