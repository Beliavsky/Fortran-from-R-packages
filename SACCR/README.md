# SACCR modern Fortran port

This project translates the computational code of the R package `SACCR` 3.4
into modern Fortran. It preserves the upstream GPL-3 license, uses an FPM
package layout, retains the original R sources and data for attribution, and
uses the previously translated `Trading` package as a local FPM dependency.

The R package's `data.tree`, JSON, package-attachment, and reflection
infrastructure is not reproduced. The calculation tree is represented by
strongly typed Fortran result records containing the same asset-class,
hedging-set, trade-level, replacement-cost, PFE, and EAD information.

## Implemented calculations

- exposure at default, projected future exposure, replacement cost, and the
  collateral multiplier;
- standard, simplified, and original-exposure-method calculation modes;
- supervisory factors, correlations, and option volatilities;
- adjusted notionals, maturity factors, and supervisory deltas;
- FX, interest-rate, credit, commodity, equity, and other-exposure add-ons;
- IRD maturity-bucket aggregation;
- credit and equity systematic/idiosyncratic aggregation;
- commodity-type and commodity-sector aggregation;
- basis and volatility hedging supersets with their 0.5 and 5 multipliers;
- CSA allocation, collateral treatment, unmargined EAD caps, and
  multi-counterparty portfolios;
- CCR methodology eligibility;
- the FX credit-protection example;
- native CSV calculation and exposure-summary export;
- all computational examples supplied by the R package.

## Dependency layout

The release is self-contained:

```text
dependencies/trading/
```

is the previously translated modern Fortran `Trading` package. The root
`fpm.toml` refers to it through a local path dependency.

## Build with FPM

```text
fpm build
fpm test
fpm run saccr_demo
```

## Build with GNU Fortran

On Unix-like systems:

```text
./build_gfortran.sh
```

On Windows:

```text
build_gfortran.bat
```

The scripts build the bundled `Trading` dependency, run the SACCR tests, and
run the demonstration program.

## Minimal example

```fortran
program example
  use saccr, only : dp, portfolio_result_t, example_ird
  implicit none

  type(portfolio_result_t) :: result

  call example_ird(result)
  write(*, '(f12.2)') result%total_ead
end program example
```

The result is approximately `569.47`, matching the package's documented Basel
IRD example.

## Main modules

- `saccr`: umbrella module
- `saccr_types`: typed calculation and breakdown results
- `saccr_supervisory`: default and CSV supervisory data
- `saccr_core`: EAD, PFE, RC, trade-level add-ons, and methodology
- `saccr_addon`: asset-class and hedging-set aggregation
- `saccr_portfolio`: CSA allocation and portfolio orchestration
- `saccr_io`: CSV calculator and summary writers
- `saccr_examples`: translated package examples

## Data

The upstream SACCR CSV files are in `data/`. The full upstream R package is
retained in `original-r/`. Programs using relative data paths should run from
the package root.

## Compatibility notes

The public R calculation tree is represented by `portfolio_result_t`,
`exposure_result_t`, `addon_result_t`, `asset_class_result_t`,
`hedging_set_result_t`, and `single_trade_addon_t`. See `docs/API_MAP.md` for
the API crosswalk and `docs/PORTING_NOTES.md` for documented corrections to
several upstream edge cases.

## License

The upstream package declares GPL-3. This translation is distributed under
`GPL-3.0-only`. The bundled Trading translation is also GPL-3.0-only.
