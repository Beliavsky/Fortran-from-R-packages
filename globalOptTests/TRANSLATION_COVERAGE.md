# Translation coverage

## Fully translated

The package contains 50 C benchmark functions and three R metadata helpers.
All are translated:

- all 50 objective functions in `src/objFun.c`
- `goTest`
- `getDefaultBounds`
- `getProblemDimen`
- `getGlobalOpt`

The R `.C` registration layer is replaced by ordinary Fortran procedure calls.
There is no plotting code in the upstream package.

## Source-level quirks

The goal is compatibility with the upstream package, not correction of every
historical benchmark formula.

### Hartman3 and Hartman6

The C source declares four rows in `a` and `p`, but loops five times.  The
coefficient array `c` is declared length five with only four initializers, so
the fifth coefficient is zero.  Accessing row five of `a` and `p` is undefined
C behavior even though the intended fifth contribution is zero.

The Fortran translation evaluates the four defined Hartmann terms.  At the
reference test points this agrees with the compiled upstream C result while
removing the out-of-bounds access.

### Gulf

The C loop starts at `j=0` and evaluates `log(0)`.  For the documented domain
the corresponding residual has limiting value zero and contributes nothing to
the objective.  The Fortran version starts at `j=1`, explicitly omitting that
zero limiting contribution and avoiding an invalid logarithm.

### PowellQ

The source formula contains

```text
(x[0] + 10*x[0])^2
```

rather than the more familiar Powell expression involving `x[1]`.  The
translation preserves the source expression exactly.

### Easom bounds

The upstream R metadata gives upper bounds `(10, 2)`, even though the usual
Easom minimizer is near `(pi, pi)`.  The translation preserves the published
metadata exactly.

### Objective pi versus metadata pi

The C objectives use the macro `3.14159265359`.  The R bounds use R's full
`pi`.  The translation keeps these two constants separate so both behaviors
are reproduced closely.

## Intentionally omitted

- R `.C` marshalling and package registration
- R argument matching/error text
- package printing/documentation infrastructure

These are not numerical algorithms.
