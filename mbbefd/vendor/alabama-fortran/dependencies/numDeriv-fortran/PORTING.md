# Porting notes

## Interface changes required by Fortran

R passes arbitrary functions and additional arguments dynamically. The Fortran
translation instead uses explicit abstract callback interfaces. Additional
parameters should be stored in a module or incorporated into a derived-type
procedure design by an application.

Complex-step routines are separate entry points (`grad_complex`,
`jacobian_complex`, and `hessian_complex`) because Fortran statically
distinguishes real and complex callback signatures.

R represents an unspecified `side` element with `NA`. Fortran uses integer
`0` for the same centered-derivative meaning; `+1` and `-1` retain their R
meanings.

The R `genD` list becomes `type(gend_result)`. Its `D` member is named `dmat`
to avoid visual ambiguity with the option `d`.

## Numerical behavior

The default algorithms, steps, and extrapolation structure follow the R code.
The extrapolation factor is generalized from the upstream hard-coded `4^m` to
`v^(2m)`, which is identical for the default `v=2` and mathematically correct
for other reduction factors.

The upstream package documents an `r=1` indexing bug. This translation handles
`r=1` as a single central-difference level without extrapolation.

Callbacks returning nonfinite values produce `nd_nonfinite_value` rather than
an R exception. Vector callbacks are checked on every evaluation to ensure
that their result length remains constant.

## Omitted infrastructure

S3 generic dispatch, R list/class construction, localization files, package
help machinery, and vignette build artifacts are language-specific and are not
translated. The computational behavior of all four exports is retained.
