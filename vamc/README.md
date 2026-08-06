# vamc-fortran

Modern Fortran 2018 translation of the computational code in the R package
`vamc` 0.2.1, a Monte Carlo valuation framework for variable annuities.

## Implemented

- Swap-curve bootstrapping with 30/360, ACT/360, ACT/365, and ACT/ACT day counts.
- General and New York business-day calendars and all five upstream rolling conventions.
- Log-linear discount-factor interpolation and extrapolation.
- Correlated multivariate Black-Scholes index scenarios and fund mapping.
- Synthetic variable-annuity portfolio generation.
- Mortality factors using uniform distribution of deaths within each attained age.
- All 19 rider families used by the upstream package:
  `DBRP`, `DBRU`, `DBSU`, `ABRP`, `ABRU`, `ABSU`, `IBRP`, `IBRU`,
  `IBSU`, `MBRP`, `MBRU`, `MBSU`, `WBRP`, `WBRU`, `WBSU`, `DBAB`,
  `DBIB`, `DBMB`, and `DBWB`.
- Single-policy and portfolio valuation over one or many scenarios.
- Historical aging of one policy or an entire portfolio.

The R data-frame and three-dimensional-array interfaces are represented by typed
Fortran derived types and ordinary allocatable arrays.

## Build

```text
fpm build
fpm test
fpm run
```

GNU Make alternatives:

```text
make check
make optimized
make MODE=optimized FLAGS='-std=f2018 -Wall -Wextra -Werror -pedantic -O3' example
```

## Minimal curve example

```fortran
use vamc

real(dp) :: rates(8)
integer :: tenors(8)
type(date_type) :: curve_date
type(yield_curve_type) :: curve
type(status_type) :: status

rates = [0.69_dp, 0.77_dp, 0.88_dp, 1.01_dp, 1.14_dp, 1.38_dp, &
         1.66_dp, 2.15_dp] * 0.01_dp
tenors = [1, 2, 3, 4, 5, 7, 10, 30]

call parse_date('2016-02-08', curve_date, status)
call build_curve_named(rates, tenors, 6, 'Thirty360', 6, 'ACT360', &
                       'NY', 'Modified_Foll', curve_date, 2, 'Thirty360', &
                       curve, status=status)
```

The resulting forward-rate vector reproduces the upstream published fixture:

```text
0.006912035  0.008520036  0.011060713  0.014146403  0.016846310
0.020451281  0.024514485  0.032098279  0.032098279
```

## Source-compatible and corrected behavior

`valuation_options_type` preserves upstream rider timing and indexing by default.
Set its fields to `.false.` to use corrected maturity timing, anniversary-benefit
renewal, maturity-benefit scalar handling, and historical-scenario indexing.

The index simulator similarly preserves the upstream Cholesky-row drift adjustment
by default. Pass `source_compatible_drift=.false.` to use half of each covariance
diagonal. Calendar functions can disable the upstream user-holiday inversion through
`source_compatible_holidays=.false.`.

See `PORTING_NOTES.md` for details.

## License

The upstream package declares GPL-2. This port is distributed as
GPL-2.0-only. The complete upstream source and original uploaded archive are
retained under `upstream/` for provenance.
