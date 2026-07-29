# Porting notes

## Source package

- Package: Risk
- Version: 1.0
- Authors: Saralees Nadarajah and Stephen Chan
- Date: 2017-06-05
- Original license: GPL (>= 2)
- Original computational source: `R/Risk.R`

The original R file contains 26 exported numerical functions and no compiled
code. All 26 are represented in this Fortran library.

## R-to-Fortran mapping

R selects a distribution dynamically from a string such as `"norm"` and then
calls `dnorm`, `pnorm`, and `qnorm`. Fortran replaces this with an abstract
`continuous_distribution` type having deferred `pdf`, `cdf`, `quantile`, and
support-bound methods.

This retains the package's open-ended distribution support without relying on
strings, global state, or external R libraries. A `callback_distribution`
provides the shortest route for user code, while normal, lognormal, uniform,
exponential, logistic, and Student t implementations are included directly.

R's adaptive `integrate` and `uniroot` calls are replaced by self-contained
adaptive Simpson integration, transformations for infinite intervals, and a
bracketed bisection solver.

## Preserved source semantics

Several names do not match the most common modern definitions. The Fortran
procedures intentionally follow the executable R source so that translated
results remain source-compatible.

### `esg`

The R code computes

```text
(1 / alpha) integral from 0 to alpha of quantile(p) dp
```

which is the lower-tail conditional mean for a continuous distribution. The
associated article describes an upper-tail expected shortfall instead. The
Fortran `esg` preserves the R implementation by evaluating the equivalent
lower partial first moment.

### `bvar`

The R code evaluates the same lower partial first moment, but with an explicit
lower integration limit `a`. This behavior is preserved.

### `expp`

The R code balances squared lower and upper partial deviations and finds the
zero of that balance. This differs from the first-order condition of the
standard asymmetric least-squares expectile objective. The source equation is
preserved. The positive-part operation is evaluated pointwise.

### `tcm`

The original code solves a CDF equation numerically. For a continuous
probability distribution it is exactly the quantile at `(1 + alpha) / 2`, so
the Fortran implementation uses that identity directly.

### `saring3`

The formula contains `aa / (2 * (aa - 1))`. Although the original example uses
`aa = 1`, the expression is undefined there. The Fortran routine returns a
quiet NaN instead of triggering a floating-point division-by-zero exception.

## Numerical details

- Infinite endpoints are represented by values whose magnitude is at least
  half of `huge(1.0_dp)`.
- Whole-line integration uses a tangent transformation.
- Semi-infinite integration uses a rational transformation.
- Default absolute and relative integration tolerances are `1e-9`.
- Expectile root calculations use tighter internal partial-moment tolerances
  and bracket expansion when necessary.
- Invalid domains return quiet NaNs rather than stopping the program.

## Omitted infrastructure

There was no plotting, data, S3/S4 class system, compiled extension, or package
vignette in the attached source. R documentation markup and dynamic dispatch
were replaced by Markdown documentation and the Fortran type system.
