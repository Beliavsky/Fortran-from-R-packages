# Porting notes

## Numerical kernel

The supplied package contains the Morales-Nocedal L-BFGS-B 3.0 Fortran kernel,
modified to use integer task codes and R-compatible printing hooks. It was
converted to free-form source and placed in `lbfgsb3_core_mod`. The task codes,
workspace layout, correction memory, Cauchy-point calculation, subspace solve,
and More-Thuente line search are retained.

For close source traceability, labeled `do` termination remains in parts of the
kernel. This is an obsolescent but valid Fortran 2018 feature; all public wrapper
code uses structured modern constructs.

## Removed dependencies

The R package obtains `dpofa` and `dtrsl` from R's bundled LINPACK interface.
This port includes compatible self-contained Cholesky and triangular-solve
implementations. No external BLAS, LAPACK, R, Rcpp, or C compiler is required.

## Bounds

A bound array may have length one or `size(x)`. IEEE infinities indicate missing
bounds internally. Equal lower and upper values fix a parameter.

## Stopping rules

The original R wrapper's `maxit` is effectively checked against objective/
gradient requests, despite its name. The Fortran control type exposes both
`max_evaluations` and `max_iterations`. The default parameter-change rule
matches the package defaults (`reltol=1e-6`, `abstol=0`). Set both to zero to
rely only on native L-BFGS-B convergence tests.

## Finite differences

`lbfgsb_minimize_fd` uses centered differences where both directions are
feasible and one-sided differences at active bounds. Its function count reports
actual objective calls rather than only L-BFGS-B gradient requests.

## Printing

The original R-specific `intpr`/`dblepr` hooks are replaced by Fortran output
routines. They are disabled for `iprint < 0`; the higher-level `trace` control is
recommended for portable progress reporting.
