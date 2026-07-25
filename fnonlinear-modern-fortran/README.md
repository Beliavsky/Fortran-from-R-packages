# fNonlinear modern Fortran

A modern Fortran translation of the computational routines in the Rmetrics
`fNonlinear` package. Plotting, R classes, formula handling, and time-series
metadata are intentionally excluded; the numerical arrays and statistics that
underlie those interfaces are available directly.

## Implemented numerical surface

### Chaotic maps and ODE integration

- Tent map
- Henon map
- Ikeda map
- Logistic map
- Lorenz/Lorentz system
- Rossler/Roessler system
- Generic fixed-step fourth-order Runge-Kutta integration over an explicit
  increasing time grid

### Nonlinear time-series statistics

- Delay embedding using regular or explicit lag sets
- Matrix delay embedding
- Average mutual information
- False-nearest-neighbor fractions and counts
- Direct, box-index, and automatic neighbor search
- Recurrence distance and threshold matrices
- Space-time separation quantiles
- Correlation integral at one radius
- Multi-radius, multi-dimension correlation-integral curves
- Nearest-neighbor Lyapunov stretching paths
- Linear Lyapunov-slope estimation

### Hypothesis tests

- Brock-Dechert-Scheinkman test using the original asymptotic variance formula
- White neural-network nonlinearity test
- Terasvirta neural-network nonlinearity test
- Runs test after discarding exact zeros
- `ts_test` method dispatcher for BDS, White, and Terasvirta tests

The test result classes used by R are replaced by plain Fortran derived types.

## Building

The validated Unix build requires GNU Fortran, LAPACK, and BLAS:

```sh
make check
make check-opt
```

`make check` enables bounds and runtime checking. `make check-opt` uses `-O2`.
Both treat compiler warnings as errors.

An `fpm.toml` manifest is included:

```sh
fpm test
fpm run demo_fnonlinear
fpm run analyze_csv -- data/logistic_sample.csv auto
```

`fpm` was not installed in the validation environment, so the manifest is
included but not claimed as tested.

## CSV application

The analyzer accepts either one numeric column or `Date,Value`/`Index,Value`
rows with an optional header:

```sh
./build/debug/analyze_csv data/logistic_sample.csv auto
./build/debug/analyze_csv data/logistic_sample.csv box
./build/debug/analyze_csv data/logistic_sample.csv direct
```

It reports mutual information, a correlation integral, false-neighbor
fractions, BDS statistics, White and Terasvirta p-values, a runs-test p-value,
and an early Lyapunov stretching slope when enough neighbors are available.

## Numerical differences from R

- The BDS implementation evaluates pairs directly instead of using the
  original bitmap acceleration. It preserves the shared effective sample,
  Dechert `k` statistic, correlation integrals, and normalization formula.
- White-test random weights use a reproducible Park-Miller generator rather
  than R's RNG. LAPACK eigendecomposition replaces `prcomp`; test statistics
  follow the same regressions and degrees-of-freedom formulas.
- Neural-test p-values use self-contained incomplete-gamma and incomplete-beta
  evaluations.
- The TISEAN-derived neighbor routines return plain arrays. The box index is a
  tested performance extension and is not represented as an R option in the
  original package.
- RK4 follows the package's fixed-step formula but does not reproduce R object
  attributes or external solver state.

Exact R random streams, PCA signs, floating-point tie behavior, and
iteration-for-iteration output are not claimed.

## Exclusions

- Plotting and graphical arguments
- S4 `fHTEST`, `timeSeries`, and `ts` objects
- R formulas, model frames, printing, and descriptions
- GUI or interactive infrastructure

No other self-contained numerical algorithms were found in the attached
`fNonlinear` source.

## License

The original package declares `GPL (>= 2)`. This translation uses
`SPDX-License-Identifier: GPL-2.0-or-later` in every Fortran source and includes
the GNU GPL version 2 text in `LICENSE`.
