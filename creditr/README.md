# creditr-fortran

This application is based on the **ISDA CDS Standard Model (version 1.0), developed and supported in collaboration with Markit Group Ltd.**

`creditr-fortran` is a modern Fortran/FPM translation of the computational algorithms in the R package `creditr` 0.6.2. It provides a self-contained CDS date engine, interest-rate curve bootstrap, one-name CDS hazard calibration, clean and dirty upfront valuation, spread/upfront conversion, and the package's risk measures.

## Implemented computational scope

- CDS trade, step-in, value, accrual, coupon, maturity, backstop, and curve-base dates
- USD, EUR, and JPY money-market and swap conventions
- Weekend, modified-following, following, and package-compatible holiday adjustments
- ACT/360, ACT/365F, and 30/360 year fractions
- Money-market and par-swap zero-curve bootstrap
- Flat-forward interpolation in logarithmic discount factors
- Quarterly CDS premium schedules
- Premium-leg accrued-on-default integration
- Protection-leg default integration
- Clean-spread flat-hazard bootstrap
- Clean and dirty upfront valuation
- Conventional spread from upfront
- Spread DV01, CS10, recovery risk 01, and interest-rate DV01
- Simple spread/default-probability/recovery/PV01 formulas
- CSV rate-curve input

The original package tree, including the ISDA C source and the bundled `rates.RData`, is retained under `original/` for provenance.

## Build with FPM

```text
fpm test
fpm run
fpm run --example spread_upfront_roundtrip
```

GNU Fortran users can run the included build script when FPM is unavailable:

```text
./scripts/build_gfortran.sh debug
./scripts/build_gfortran.sh release
```

On Windows:

```text
scripts\build_gfortran.bat debug
scripts\build_gfortran.bat release
```

## Minimal example

```fortran
use creditr

type(rate_quote_t), allocatable :: quotes(:)
type(conventions_t) :: conventions
type(zero_curve_t) :: curve
type(cds_contract_t) :: contract
type(cds_result_t) :: result
integer :: status

call read_rate_quotes_csv('data/usd_2014_04_15.csv', quotes, status)
call add_conventions('USD', conventions, status)
call build_zero_curve(make_date(2014, 4, 17), quotes, conventions, curve, status)

contract%trade_date = make_date(2014, 4, 15)
contract%maturity = make_date(2019, 6, 20)
contract%use_maturity = .true.
contract%spread_bps = 243.28_dp
contract%coupon_bps = 100.0_dp
contract%currency = 'USD'

result = price_cds(contract, curve, status, quotes)
```

Passing `quotes` to `price_cds` makes interest-rate DV01 rebuild the curve after adding one basis point to every market quote, matching `creditr::IR_DV01`. Without the optional quotes, the routine applies a parallel continuous-zero-rate bump to the supplied curve.

## Dependencies

The library has no external numerical dependency. It uses only standard Fortran. The original R package's `RCurl`, `XML`, `zoo`, `xts`, Rcpp, and online Markit download machinery are not required.

## Numerical validation

The tests include:

- exact package date and convention cases;
- all 64 nodes of the representative USD curve structure, with fixed ISDA discount-factor references;
- four one-name CDS valuations spanning 243 to 12,355 basis points;
- hazard-rate, spread-DV01, and recovery-risk references generated from the retained ISDA C engine;
- clean/dirty accrual identities;
- spread/upfront round trips; and
- quote-bumped interest-rate sensitivity.

See `TESTING.md` and `REFERENCE_GENERATION.md`.

## Scope differences

The numerical valuation code is translated. R S4 classes, formatted printing, HTML/XML downloads, and automatic lookup in `rates.RData` are not reproduced. Rate quotes are supplied explicitly or loaded from CSV. More details are in `PORTING.md`.

## License

The original R work is GPL-3.0-only. Algorithms derived from the ISDA source are also subject to the **ISDA CDS Standard Model Public License 1.0**. The complete combined original license is in `LICENSE`; separate copies are in `LICENSES/`. See `NOTICE.md` before redistribution.
