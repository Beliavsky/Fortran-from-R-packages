# Translation status

Upstream package: `quadform` 0.0-4.

The package contains no compiled code and no non-base computational dependency. Its full computational R source is in `R/quadform.R`; all exported numerical routines in that file have Fortran counterparts in this translation.

No plotting code exists in the package.

## Deliberate Fortran differences

- R matrices may sometimes be created implicitly from vectors; the core Fortran API is explicitly matrix-oriented.
- R's `solve()` signals an error for singular systems. The Fortran inverse-form routines instead return an empty matrix and an optional nonzero `info` code, which is friendlier for a numerical library.
- R object attributes, dimnames, recycling, and other language-level semantics are not reproduced.
