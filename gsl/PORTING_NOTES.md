# Porting notes

## Architecture

The upstream R package delegates essentially all numerical special-function work to GNU GSL. A literal independent Fortran rewrite of those algorithms would no longer be a translation of this package; it would be a reimplementation of GNU GSL itself. The Fortran port therefore keeps GNU GSL as the numerical backend and translates the R-facing interface into modern Fortran.

The implementation has three layers:

1. `gsl_special.f90`: typed Fortran wrappers for the package's 239 non-R compiled special-function entry points.
2. `src/c_shim/*.c`: the upstream thin GSL shims, excluding R registration and R external-pointer code.
3. Native Fortran object wrappers for RNG/QRNG state and R-only numerical conveniences.

## Case sensitivity

Fortran identifiers are case-insensitive. The R package uses case to distinguish cylindrical and spherical Bessel functions. The port uses explicit `cyl`, `sph`, and `mod` components in those names.

## Vector recycling

R's `process.args()` recycles shorter arguments to the longest length and preserves attributes. The Fortran API deliberately does not imitate this language-level behavior. Input arrays should already have compatible lengths. This makes shape mismatches explicit and avoids hidden allocations.

## Error semantics

The low-level wrappers preserve GSL status and estimated-error outputs. `strictify()` implements the R package's default idea of replacing unsuccessful values by NaN, but it is explicit rather than automatic.

## Legendre raw/e variants

The upstream C layer contains both direct-return and `_e` variants for several Legendre functions. The R API uses the `_e` versions so that status/error estimates are available. In Fortran, those are the primary names; direct-return wrappers are retained with a `_raw` suffix.

## Polynomial evaluation

`gsl_poly` is implemented directly in Fortran with Horner's rule. The upstream C-wrapper version remains reachable as `gsl_poly_c` for interface parity.

## RNG ownership

R external pointers provide automatic finalizers. Plain Fortran derived-type assignment can duplicate opaque C pointers, so automatic finalization would risk double frees. The port therefore requires explicit `rng_free` / `qrng_free`. Cloning creates an independent GSL state.

## multimin

The package exports `multimin`, `multimin.init`, `multimin.iterate`, and related names, but the current R implementation immediately stops and states that the functions are temporarily removed. No active numerical behavior exists to translate in version 2.1-9, so these are documented rather than recreated.
