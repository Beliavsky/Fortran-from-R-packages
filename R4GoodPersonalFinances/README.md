# R4GoodPersonalFinances-fortran

A modern Fortran/FPM translation of the reusable computational core of the R
package `R4GoodPersonalFinances` 1.2.0.9000.

The library implements mortality and longevity calculations, retirement ruin,
personal-finance utility functions, tax-aware portfolio construction,
household timelines, correlated return generation, and lifetime household
balance-sheet simulation. It is self-contained and has no external numerical
dependencies.

## Build

```text
fpm build
fpm test
fpm run
```

GNU Fortran scripts are also supplied:

```text
scripts/build_checked.sh
scripts/build_optimized.sh
```

On Windows with `gfortran` on `PATH`:

```bat
scripts\build_checked.bat
scripts\build_optimized.bat
```

## Main modules

- `r4gpf_mortality`: Gompertz survival, calibration, joint-life fitting, life
  expectancy, incomplete gamma functions, and retirement ruin.
- `r4gpf_finance`: purchasing power, present values, CRRA utility, certainty
  equivalents, Merton allocation, and effective tax rates.
- `r4gpf_portfolio`: portfolio records, moments, tax-aware and total-net-worth
  optimization, and correlated Gaussian returns.
- `r4gpf_household`: dates, household members, events, joint survival,
  timelines, and callback-based cash-flow streams.
- `r4gpf_simulation`: optimal discretionary spending and deterministic or
  Monte Carlo lifetime simulations.
- `r4good_personal_finances`: umbrella module exporting the public API.

Fortran arrays are explicit and strongly typed. Asset vectors have length
`n_assets`; return matrices are `periods x assets`; allocation matrices are
`assets x periods`.

## Examples

Five focused examples are under `example/`. The default executable
`demo_r4good_personal_finances` constructs a household and portfolio, computes
longevity and ruin measures, and runs a short lifetime simulation.

## Scope

The numerical functions are direct translations where practical. Dynamic R
formula evaluation is replaced by procedure callbacks, and lifetime simulation
accepts explicit income and spending vectors. See `API_MAP.md` and
`PORTING_NOTES.md` for exact coverage and adaptations.

The Shiny application, plotting, report rendering, cache/parallel orchestration,
R6/S3 infrastructure, embedded `.rda` life tables, and package fonts/images are
not part of the Fortran runtime.

## License

MIT. See `LICENSE.md`, `NOTICE.md`, and the retained upstream snapshot.
