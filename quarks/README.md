# quarks-fortran

Modern Fortran 2018 translation of the computational core of the R package
**quarks 1.1.6**, packaged for the Fortran Package Manager (FPM).

The library calculates and backtests Value at Risk (VaR) and Expected Shortfall
(ES) using plain historical simulation, age weighting, volatility weighting,
and filtered historical simulation. The GARCH path uses the included modern
Fortran translation of `rugarch`; EWMA paths are dependency-free.

## Implemented routines

- `ewma`: exponentially weighted conditional variance.
- `hs`: plain or age-weighted historical VaR and ES.
- `vwhs`: EWMA- or GARCH-weighted historical simulation.
- `fhs`: bootstrap filtered historical simulation.
- `rollcast`: moving-window one-step VaR and ES forecasts.
- `cvgtest`: Kupiec/Christoffersen coverage and independence tests.
- `trftest`: Basel traffic-light cumulative breach probability.
- `lossfun`: the four upstream ES loss functions.
- `plop` and `plop_time_varying`: portfolio return/loss aggregation.

## Build

```console
fpm build
fpm test
fpm run
```

The project contains a relative FPM path dependency at
`vendor/rugarch-modern-fortran`, so no system-installed package is required.

GNU Fortran scripts are also supplied:

```console
./scripts/run_checked.sh
./scripts/run_optimized.sh
```

On Windows:

```bat
scripts\run_checked.bat
scripts\run_optimized.bat
```

## Basic example

```fortran
use quarks
real(dp) :: returns(1000)
type(risk_result) :: risk

! Fill returns first.
risk = hs(returns, p=0.975_dp, method=method_age, lambda=0.98_dp)
print *, risk%var, risk%es
```

## Data orientation

Return series are one-dimensional arrays ordered from oldest to newest.
Portfolio matrices are `observations x assets`. Rolling forecasts use only
observations available before each forecast date.

## Important adaptations

The experimental `smoots`-based unconditional-scale options in the R package
are represented by a native Gaussian-kernel local-scale estimator. The default
`smooth_none` path is a direct translation.

The R package forwards arbitrary `rugarch::ugarchspec` arguments. The compact
Fortran `vwhs` and `fhs` interfaces expose the translated Gaussian sGARCH(1,1)
path and EWMA. Advanced users can call the vendored `rugarch` API directly.

See `API_MAP.md` and `PORTING_NOTES.md` for details.

## License and status

This combined project is distributed under GPL-3.0-only. It is an independent,
experimental translation and is not an official release by the original
package author. Verify risk estimates independently before production use.
