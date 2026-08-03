# fincal-fortran

Modern Fortran translation of the computational routines in **FinCal 0.6.3**, an R package for time-value-of-money calculations, yield conversions, financial ratios, portfolio statistics, depreciation, and inventory costing.

The translation is self-contained and uses the Fortran Package Manager (FPM). Plotting routines and network-dependent historical-price download functions are intentionally omitted.

## Build and test

```console
fpm build
fpm test
fpm run demo_fincal
```

A compiler-only test script is also supplied:

```console
./run_gfortran_tests.sh
```

On Windows:

```console
run_gfortran_tests.bat
```

## Use

```fortran
use fincal
```

All calculations use `real(dp)`, where `dp = kind(1.0d0)`.

```fortran
real(dp) :: value, rate
integer :: status

value = fv(0.07_dp, 10.0_dp, present_value = 1000.0_dp, payment = 10.0_dp)
rate = irr([-5.0_dp, 1.6_dp, 2.4_dp, 2.8_dp], status)
```

See [API.md](API.md) for the complete mapping from R names to Fortran names.

## Scope

Translated computational functions:

- time value of money: `fv`, `pv`, annuities, uneven cash flows, payments, periods, discount rates, NPV, IRR, and perpetuities
- interest and money-market conversions: EIR, EAR, HPR, BDY, MMY, BEY, and continuous/nominal rates
- portfolio/statistical functions: geometric and harmonic means, weighted return, TWRR, Sharpe ratio, Roy safety-first ratio, coefficient of variation, and sampling error
- accounting functions: EPS, diluted EPS, weighted shares, treasury-stock-method issuable shares, straight-line and double-declining depreciation, and FIFO/LIFO/WAC inventory costing
- liquidity, solvency, and profitability ratios

Omitted functions:

- `candlestickChart`, `lineChart`, `lineChartMult`, and `volumeChart` because they are plotting-only
- `get.ohlc.google`, `get.ohlc.yahoo`, `get.ohlcs.google`, and `get.ohlcs.yahoo` because they are network/data-retrieval wrappers rather than computational algorithms

## Numerical notes

`irr` and `discount_rate` use a self-contained transformed-rate scan followed by bisection. `irr2` retains the original package's thresholded linear search, including support for negative rates above -100 percent. See [PORTING.md](PORTING.md) for details.

## License

The original FinCal package declares `GPL (>= 2)`. This translation is therefore distributed under **GPL-2.0-or-later**. Original computational R sources and package metadata are retained under `original/` for attribution and provenance.
