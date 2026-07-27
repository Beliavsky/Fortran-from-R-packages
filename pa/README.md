# pa-modern-fortran

A modern Fortran translation of the non-plotting computational core of the R package `pa`, which performs equity-portfolio performance attribution using Brinson and regression-based methods.

The project uses plain arrays and derived types rather than R data frames and S4 classes. It is licensed under GPL-2.0-only, matching the original package's `GPL-2` declaration.

## Implemented numerical functionality

### Brinson attribution

- Single-period portfolio and benchmark weights by category.
- Weight-adjusted portfolio and benchmark category returns.
- The four Brinson quadrant returns `q4`, `q3`, `q2`, and `q1`.
- Category allocation, selection, and interaction effects.
- Aggregate allocation, selection, interaction, and active return.
- Multi-period attribution with a common category panel.
- Arithmetic aggregation.
- Geometric aggregation.
- Linking-coefficient aggregation using the original package formulas.

### Regression attribution

- No-intercept least-squares factor-return estimation using LAPACK `DGELSY`.
- Active factor exposures from portfolio-minus-benchmark weights.
- Factor contributions and residual contribution.
- Portfolio, benchmark, and active returns.
- Single- and multi-period results.
- Arithmetic, geometric, and linking aggregation.
- Conceptual-variable contribution aggregation through explicit column ranges.

### Design matrices and exposures

- Mixed categorical and numeric design-matrix construction.
- R-compatible no-intercept treatment coding: every level of the first categorical predictor, and all but the first level of later categorical predictors.
- R-style average ranks for ties.
- Quintile assignment using `ceiling(5 * rank / n)`, matching the original exposure methods.
- Categorical and continuous exposure calculations.
- Single- and multi-period portfolio, benchmark, and active exposures.

## Build

GNU Fortran, LAPACK, and BLAS are required.

```sh
make debug
make release
```

The debug build uses full runtime checking and the release build uses optimization. Both treat compiler warnings as errors.

An `fpm.toml` file is included:

```sh
fpm test
fpm run demo_pa
```

`fpm` was not available in the validation environment, so the Makefile builds are the tested build route.

## CSV application

The input columns are:

```text
period,category,benchmark,portfolio,return,value,growth
```

Examples:

```sh
fit_csv data/example_attribution.csv brinson geometric
fit_csv data/example_attribution.csv brinson linking
fit_csv data/example_attribution.csv regression arithmetic
```

The regression CSV mode uses `category` as a categorical factor and `value` and `growth` as numeric factors.

## Numerical differences and safeguards

- R formula and factor processing is replaced by an explicit array design-matrix builder.
- LAPACK rank-revealing least squares replaces R's `lm` implementation. Coefficients agree for full-rank designs, but rank-deficient coefficient choices can differ while fitted values remain equivalent.
- Integer period and category identifiers replace R date/factor metadata.
- Zero category weights produce a category return of zero rather than propagating an undefined division.
- Degenerate linking cases are guarded and reported through the result status field.

## Excluded R infrastructure

- S4 classes and methods.
- Data-frame and formula evaluation.
- Plotting and `ggplot2`/`grid` graphics.
- Console `show` and formatted `summary` methods.
- R factor labels, date classes, row names, and column names.
- Packaged R `.RData` datasets.

The numerical arrays underlying the excluded displays are available through the Fortran result types.
