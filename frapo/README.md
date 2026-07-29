# FRAPO modern Fortran port

This project translates the computational algorithms of **FRAPO 0.4-2** to
modern Fortran and packages them as an FPM library.

The port includes return and trend transformations, portfolio risk measures,
tail-dependence estimation, matrix square roots, long-only portfolio
construction, equal-risk-contribution optimization, and the four drawdown
portfolio formulations. It has no R dependency.

## Implemented algorithms

- `capser`: lower and upper capping of vectors or matrices
- `returnseries`: continuous or discrete returns, optional compounding,
  percentage scaling, and trimming
- `returnconvert`: continuous/discrete return conversion
- `trdbinary`, `trdbilson`, `trdes`, `trdhp`, `trdsma`, and `trdwma`
- `mrc`, `dr`, `cr`, and `rhow`
- `tdc`: empirical and EVT lower- or upper-tail dependence coefficients
- `sqrm`: real symmetric matrix square root
- `PGMV`, `PMD`, `PMTD`, and `PERC`
- `PMaxDD`, `PAveDD`, `PCDaR`, and `PMinCDaR`

Both descriptive Fortran names and the original short FRAPO names are public.
Fortran is case-insensitive, so `pgmv` is the counterpart of R's `PGMV`.

## Dependencies

- A Fortran 2018 compiler
- BLAS
- LAPACK
- FPM, optionally

No external optimization package is needed. The R dependencies `cccp` and
`Rglpk` are replaced by self-contained numerical solvers.

## Building with FPM

```text
fpm build
fpm test
fpm run
```

The manifest links against `lapack` and `blas`.

## Reproducible GNU Fortran build

```text
./run_tests.sh strict
./run_tests.sh optimized
```

The script compiles the library, runs every test, and compiles every application
and example.

## Minimal example

```fortran
use frapo
implicit none

real(dp) :: returns(100, 4)
type(portfolio_result) :: solution

! Fill returns(:, :) here.
solution = pgmv(returns, percentage=.false.)
print *, solution%weights
```

See `API.md`, `PORTING.md`, and the programs in `example/` for details.

## Licensing and provenance

The Fortran port is distributed under **GPL-3.0-or-later**, matching the
attached FRAPO package. The complete original package tree is retained under
`original/FRAPO-master/` for provenance.

Some retained index-tracking datasets originate from J. E. Beasley's
OR-Library and are covered by the MIT-style terms reproduced in
`BEASLEY-LICENSE.txt` and in the original package.
