# fBonds Modern Fortran

A modern Fortran translation of the computational routines in the R package
`fBonds` version 3042.78.

The original package contains two numerical term-structure models:

- Nelson-Siegel
- Nelson-Siegel-Svensson (named `Svensson` in the R package)

This project translates both models, their starting-value searches, parameter
estimation, fitted curves, residuals, and fit diagnostics. Plotting and R object
infrastructure are intentionally omitted.

## Features

- Numerically stable Nelson-Siegel curve evaluation, including maturity zero
- Svensson curve evaluation using the exact forward-rate basis in the R source
- Log-spaced global starting grids for positive decay constants
- SVD least-squares starting coefficients through LAPACK `DGELSS`
- Positive decay constants enforced through log parameterization
- Nelson-Siegel sum-of-squared-errors fitting
- Svensson absolute-error fitting, matching the effective R source behavior
- Optional sum-of-squared-errors Svensson fitting
- Fitted values, residuals, objective value, SSE, MAE, and RMSE
- Convergence status, iterations, and objective-evaluation counts
- CSV command-line application
- Synthetic and original-package example tests

## Public API

```fortran
use fbonds_term_structure

type(term_structure_fit) :: fit

call fit_nelson_siegel(rate, maturity, fit)
call fit_svensson(rate, maturity, fit)
call fit_svensson(rate, maturity, fit, objective="sse")

curve = nelson_siegel_curve(maturity, parameters)
curve = svensson_curve(maturity, parameters)
```

Nelson-Siegel parameters are ordered as:

```text
beta0, beta1, beta2, tau1
```

Svensson parameters are ordered as:

```text
beta0, beta1, beta2, beta3, tau1, tau2
```

## CSV application

The input must contain two numeric columns in this order:

```text
maturity,rate
```

A header is optional.

```text
fit_csv data/example_yield.csv ns
fit_csv data/example_yield.csv svensson l1
fit_csv data/example_yield.csv svensson sse
```

## Build

GNU Fortran, LAPACK, and BLAS are required.

```text
make debug
make release
```

The debug build enables bounds and runtime checks. Both modes treat compiler
warnings as errors.

On Windows with GNU Fortran and LAPACK/BLAS available:

```text
build.bat
build.bat release
```

An `fpm.toml` manifest is included. `fpm` was not installed in the validation
environment and is not claimed as tested.

## Source-level behavior preserved or corrected

The R `Svensson` objective contains an SSE expression followed by an absolute
error expression. R returns the final expression, so the effective objective is
L1. The Fortran default preserves this behavior; `objective="sse"` is available
explicitly.

The R Nelson-Siegel grid start solves coefficients using the standard
Nelson-Siegel yield loadings but evaluates the temporary grid curve using a
different exponential basis. The Fortran grid evaluates the same standard basis
used to solve the coefficients and by the final Nelson-Siegel curve. This is a
source-level consistency correction, documented in `ORIGIN.md`.

## Exclusions and numerical differences

- Plotting, contour, image, and perspective surfaces
- R list objects and `nlminb` control structures
- `timeDate`, `timeSeries`, and `fBasics` object integration
- Exact R optimizer endpoints and iteration histories

The Fortran optimizer is bounded Nelder-Mead with log-transformed positive decay
constants. Exact optimizer equivalence is not claimed.

## License

The original package declares `GPL (>= 2)`. This translation is distributed as
`GPL-2.0-or-later`. Every Fortran source file contains an SPDX identifier and an
explicit GPL version 2-or-later notice. See `LICENSE` and `ORIGIN.md`.
