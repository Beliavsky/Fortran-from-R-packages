# Porting notes

## Array conventions

Observations are rows and series are columns. Lag coefficient arrays use
`matrix(row,column,lag)`. Forecast covariance arrays use
`covariance(row,column,horizon)`. Some routines intentionally allocate arrays
with a zero lower bound for lag zero, for example `psi(:,:,0)`.

## VARMA sign convention

The implemented model is

```text
x_t = c + sum(Phi_j x_{t-j}) + e_t - sum(Theta_j e_{t-j}).
```

This matches the residual recursion used by the upstream package. Consequently,
`psi(:,:,1) = phi(:,:,1) - theta(:,:,1)` for a VARMA(1,1).

## Estimation

VAR, VARX, regression, factor, and known-beta VECM estimators use linear least
squares. VARMA, EWMA, DCC, and BEKK estimators use local BFGS with numerical
derivatives and bounds. These routines are self-contained and prioritize a
portable reference implementation over specialized BLAS/LAPACK performance.

## Numerical linear algebra

The package includes pivoted Gauss-Jordan inversion, Cholesky factorization,
and a Jacobi symmetric eigensolver. These are suitable for examples and modest
systems. Large production models should eventually offer an optional BLAS and
LAPACK backend while retaining the current fallback.

## Attached dependency translations

The supplied fGarch, fBasics, and mvtnorm Fortran projects were inspected for
package organization, naming, and numerical conventions. They remain separate
works under their own licenses. No source from those GPL projects is included
or linked here; required distributions, random generation, and linear algebra
are independently provided in this project.

## Scope differences

R-specific date/time-series attributes, data frames, formulas, interactive
selection, console formatting, S3/list objects, plotting, and `.rda` loading are
not meaningful in a numerical Fortran library and were omitted. See
`API_MAP.md` for specialized model-specification systems that are represented
by primitives or approximations rather than exact high-level wrappers.
