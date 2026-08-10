# bvls-fortran

Modern Fortran/FPM translation of the computational core of the R package
`bvls` 1.4 and the Stark-Parker bounded-variable least-squares algorithm.

It solves

```text
min || A*x - b ||_2    subject to lower <= x <= upper.
```

The translated solver preserves the active/bound working-set algorithm,
Householder QR solve, warm-start `istate` convention, short-cycle guard,
rank-dependence rejection, feasibility interpolation, and `3*n` major-loop
limit of the original Fortran 77 implementation.

## Build

```text
fpm build
fpm test
```

## Basic use

```fortran
use bvls, only : dp, bvls_result, bvls_fit

type(bvls_result) :: fit
call bvls_fit(a, b, lower, upper, fit)
print *, fit%x
print *, fit%deviance
```

For a warm start, pass `key=1` and a previously returned `fit%istate`.

## Result fields

`bvls_result` contains the coefficients, fitted values, residuals, negative
least-squares gradient `A^T(b-Ax)`, residual sum of squares, residual norm,
`istate`, major-loop iteration count, and a status code.

## License

The R package is GPL (>= 2), and the Stark-Parker source distributed with it
was authorized for GPL version 2 or newer distribution. The original package
is retained under `original/bvls-master/`.
