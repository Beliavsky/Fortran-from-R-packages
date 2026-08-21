# Porting notes

## Structural translation

The CRAN package already contains Powell's original fixed-form Fortran
implementations. This port:

1. converts all numerical source to free-form `.f90`;
2. places the routines in `minqa_module`, giving every internal call an
   explicit module interface;
3. uses `real64` through the public `dp` kind;
4. allocates algorithm workspaces inside the public wrappers;
5. replaces the Rcpp callback, R environment, and `.Call` registration with
   a native Fortran procedure callback;
6. replaces R lists and S3 classes with `minqa_control_t` and
   `minqa_result_t`;
7. maps the original raw error codes to the public status values used by the
   R package.

The internal Powell routines retain their original variable names, branch
structure, arithmetic order, and exact-real comparisons to minimize
numerical drift from the reference implementation.

## Corrected RESCUE trial-vector evaluation

In `src/rescue.f` from the supplied package, a new trial point is assembled
in `W(1:N)`, but the subsequent objective call is:

```fortran
F = CALFUN(N, X, IPRINT)
```

`X` is neither an argument nor an array declared in `RESCUE`; under the
original implicit rules it becomes a scalar. Passing it as an `N`-element
vector is invalid and can read unrelated memory if the rescue path is
entered.

The port calls:

```fortran
f = calfun(n, w, iprint)
```

which is consistent with the immediately preceding construction of
`W(1:N)` and with Powell BOBYQA implementations.

## Defaults

The R package computes default `rhobeg` as

```text
min(0.95, 0.2 * max(abs(x)))
```

which becomes zero for an all-zero initial vector and then fails validation.
The Fortran API uses `0.5` when that computed default is zero. Explicit
positive `control%rhobeg` values are preserved exactly.

BOBYQA omitted bounds are represented internally by large finite limits.
Explicit IEEE `-Inf` lower bounds and `+Inf` upper bounds are accepted and
converted to the same representation.

## Callback limitation

The low-level algorithms use a module-level active procedure pointer to keep
the translated arithmetic kernels unchanged. Consequently, nested,
recursive, or concurrent optimizer calls in the same process are not
supported. The optimizer is otherwise self-contained and has no R or C++
runtime dependency.
