# yrnd-fortran

Modern Fortran 2018 translation of the computational code in the R package
`yrnd` 0.1.5, packaged with FPM.

The library calibrates two- or three-component lognormal mixtures to option
prices on fixed-income futures and derives risk-neutral distributions for:

- bond-futures prices;
- short-term interest-rate futures prices and rates;
- cheapest-to-deliver bond yields;
- CTD probabilities across a delivery basket; and
- correlated government-bond yield spreads.

The translated `tvm` package is bundled as a local FPM dependency and is used
for irregular cash-flow yield calculations.

## Build

```sh
fpm build
fpm test
fpm run
```

A standalone GNU Fortran validation script is also supplied:

```sh
./scripts/validate.sh
```

## Minimal example

```fortran
program example
   use yrnd
   implicit none

   real(dp) :: strikes(5), calls(5), puts(5)
   type(lognormal_mixture_t) :: source
   type(density_result_t) :: fitted
   type(date_t) :: today, option_date, futures_date
   integer :: i

   strikes = [92.0_dp, 96.0_dp, 100.0_dp, 104.0_dp, 108.0_dp]
   source%n_components = 2
   source%meanlog(1:2) = [log(97.0_dp), log(103.0_dp)]
   source%sdlog(1:2) = [0.08_dp, 0.16_dp]
   source%weight(1:2) = [0.45_dp, 0.55_dp]

   call mixture_option_prices(source, strikes, strikes, 0.03_dp, 0.5_dp, &
      option_european, calls, puts)

   today = date_t(2026, 8, 5)
   option_date = date_t(2027, 2, 5)
   futures_date = date_t(2027, 3, 15)

   call bond_future_price(calls, strikes, puts, strikes, 2, 0.03_dp, &
      dc_act_365, option_european, source%mean(), futures_date, &
      option_date, today, fitted, grid_step=0.05_dp)

   write(*, '(a,f10.5)') "Fitted mean: ", fitted%moments(1)
   write(*, '(a,f10.5)') "Fitted standard deviation: ", fitted%moments(2)
end program example
```

## Public modules

- `yrnd`: umbrella module for the complete public API.
- `yrnd_mixture`: option pricing, mixture calibration, PDF/CDF/quantiles.
- `yrnd_bonds`: coupon schedules, yields, net basis, and CTD selection.
- `yrnd_transforms`: STIR-rate, bond-yield, and spread distributions.
- `yrnd_dates`: lightweight Gregorian dates and day-count conventions.

## Scope

Plotting, tibbles/data frames, S3 metadata, Bloomberg API calls, and other
R-specific interface code are omitted. See `docs/API_MAP.md` and
`docs/PORTING_NOTES.md` for precise coverage and design differences.
