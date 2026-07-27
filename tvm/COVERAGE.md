# Computational coverage

This file maps the R package `tvm` 0.5.2 to the modern Fortran implementation.

## Exported R interfaces

| R interface | Fortran interface | Status |
|---|---|---|
| `adjust_disc` | `adjust_disc` | Implemented |
| `cft` | `cft` | Implemented |
| `npv` | `npv` | Implemented |
| `xnpv` | generic `xnpv` / `xnpv_tau` / `xnpv_dates` | Implemented |
| `irr` | `irr` | Implemented |
| `xirr` | generic `xirr` / `xirr_tau` / `xirr_dates` | Implemented |
| `pmt` | `pmt` | Implemented |
| `rate` | generic `rate` / `loan_rate` | Implemented |
| `loan` | generic `loan` / `make_loan` returning `loan_t` | Implemented |
| `cashflow` | generic `cashflow` and `loan_t%cashflows` | Implemented |
| `disc_cf` | `disc_cf` | Implemented |
| `rem` | `rem` | Implemented |
| `find_rate` | `find_rate` | Implemented |
| `rate_curve` | four typed curve constructors | Implemented |
| `[.rate_curve` | `rate_curve_t%rates` and `%rate_grid` | Implemented |
| `disc_value` | `disc_value` and `rate_curve_t%present_value` | Implemented |
| `plot.rate_curve` | numerical curve values returned to caller | Plotting excluded |

## Loan cashflow methods

| R method | Fortran behavior |
|---|---|
| `cashflow.bullet` | `loan(..., "bullet")` |
| `cashflow.french` | `loan(..., "french")` |
| `cashflow.german` | `loan(..., "german")` |
| interest grace | `grace_int` argument |
| amortization grace | `grace_amort` argument |

## Internal curve transformations

All computational transformations in `R/CurveFuncs.R` are public Fortran
procedures:

- `fut_to_zero_eff`
- `disc_to_swap` and `swap_to_disc`
- `disc_to_zero_eff` and `zero_eff_to_disc`
- `disc_to_zero_nom` and `zero_nom_to_disc`
- `disc_to_fut` and `fut_to_disc`
- `disc_to_german`
- `disc_to_french`
- `disc_to_zero_cont` and `zero_cont_to_disc`
- `eff_to_dir` and `dir_to_eff`
- nominal and effective scaling/unscaling

The R package's `get_rate_fun` functionality is represented by the type-bound
`rate_curve_t%rates` method.

## Supporting numerical code

The Fortran project adds self-contained support required to replace R runtime
services:

- bracketed bisection root solver;
- interval extension for IRR calculations;
- Fritsch-Carlson monotone piecewise-cubic Hermite interpolation;
- scalar and vector curve evaluation;
- typed result storage and validation.

## Excluded infrastructure

The following items are not computational algorithms and are not compiled:

- `ggplot2`, `reshape2`, and `scales` plotting/presentation code;
- R S3 dispatch and list representation;
- roxygen and `.Rd` documentation machinery;
- vignette rendering and `testthat` infrastructure.

The complete original source, including those files, is retained unchanged
under `original/tvm-0.5.2`.
