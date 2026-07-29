# etrm-fortran

A modern Fortran 2018/FPM translation of the computational routines in the R
package **etrm 1.0.2**, *Energy Trading and Risk Management*, by Anders D.
Sleire.

The project preserves the original MIT license and copyright attribution. It is
an independent language port and is not an official release of the R package.

## Implemented scope

The library implements all six exported computational functions from the R
package:

- `cppi`: constant-proportion portfolio insurance;
- `dppi`: dynamic-proportion portfolio insurance;
- `obpi`: option-based portfolio insurance with Black-76 delta hedging;
- `shpi`: step-hedge portfolio insurance;
- `slpi`: stop-loss portfolio insurance; and
- `msfc`: maximum-smoothness forward curves for flow-delivery contracts.

The strategy routines return a typed `strategy_result` containing market price,
trades, exposed volume, hedge position, hedge fraction, target price, and
portfolio price. `summarize_strategy` provides the original summary concepts,
including churn and first/max/min/last statistics.

The MSFC implementation constructs the quartic piecewise-polynomial constrained
system from the R source and solves its KKT equations with LAPACK. It returns
daily curve values, knots, polynomial coefficients, selected contract metadata,
and continuous-time repriced contract values.

## Deliberately omitted

- S4 classes and generic `show`, `summary`, and `plot` dispatch;
- `ggplot2` and `reshape2` plotting;
- bundled `.rda` example datasets; and
- R data-frame and date-container infrastructure.

Integer day offsets relative to the curve trade date replace R `Date` vectors.
The `day_offset` convenience function converts civil dates into those offsets.

## Build and test

The project uses the standard FPM layout and requires BLAS/LAPACK:

```text
fpm build
fpm test
fpm run
fpm run --example msfc_example
fpm run --example strategies_example
```

On systems where BLAS/LAPACK are not discovered automatically, make the
libraries visible to the system linker. The manifest declares both libraries in
its `link` list.

## Minimal strategy example

```fortran
program example_cppi
   use etrm
   implicit none

   real(dp) :: futures(6)
   type(strategy_result) :: result
   integer :: status
   character(len=:), allocatable :: message

   futures = [100.0_dp, 102.0_dp, 105.0_dp, 103.0_dp, 108.0_dp, 112.0_dp]
   call cppi(10.0_dp, futures, 0.10_dp, 0.12_dp, result, status, message)
   if (status /= etrm_ok) error stop message

   print *, result%position
   print *, result%portfolio
end program example_cppi
```

## Minimal MSFC example

```fortran
logical :: include(3)
character(len=5) :: names(3)
integer :: starts(3), ends(3), status
real(dp) :: prices(3)
type(msfc_result) :: curve
character(len=:), allocatable :: message

include = .true.
names = [character(len=5) :: "JUL", "AUG", "SEP"]
starts = [14, 45, 76]
ends = [43, 75, 105]
prices = [32.55_dp, 32.50_dp, 32.50_dp]

call msfc(include, names, starts, ends, prices, curve, status, message)
```

See `API.md`, `PORTING.md`, and `TESTING.md` for details.

## Citation and provenance

The original package requests citation of:

> Sleire, Anders D. (2022). "etrm: Energy Trading and Risk Management in R."
> The R Journal 14(1), 320-341. DOI: 10.32614/RJ-2022-013.

Original R source files and package metadata are retained under `original/` for
license and implementation provenance.
