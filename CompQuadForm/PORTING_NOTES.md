# Porting notes

## Scope

The computational surface of CompQuadForm 1.4.4 is small and concentrated in
four routines. All four are translated:

- `davies`: Davies AS 155 (`src/qfc.cpp` + `R/davies.R`)
- `farebrother`: Farebrother/Ruben AS 204 (`src/ruben.cpp` + `R/farebrother.R`)
- `imhof`: Imhof inversion integral (`src/imhof.cpp` + `R/imhof.R`)
- `liu`: Liu-Tang-Zhang approximation (`R/liu.R`)

`R/zzz.R` and `src/registerDynamicSymbol.c` are R dynamic-loading and
registration glue and are intentionally omitted.

## r_mod reuse

The supplied MIT-licensed `r_mod.f90` is reused rather than duplicating
R-compatible helpers. In particular:

- Farebrother uses `normal_cdf`.
- Imhof uses `integrate` and `integrate_result_t`.
- Liu's chi-square calculation uses `pgamma` as the central chi-square
  building block.

One package-local helper, `precise_pchisq`, was added because the supplied
`r_mod` implementation of `pchisq` intentionally uses a Wilson-Hilferty
approximation for central chi-square probabilities. That approximation causes
a visible loss of parity in Liu's method (including a case that should reduce
exactly to a single noncentral chi-square). `precise_pchisq` therefore uses
`r_mod`'s existing exact `pgamma` helper plus the standard Poisson mixture for
noncentral chi-square probabilities. No gamma/incomplete-gamma routine is
reimplemented.

## Numerical translation details

Davies is a close structural translation of the upstream AS 155 C++ code.
The original file-wide static state and `setjmp` escape are replaced by local
state captured by internal procedures and an explicit limit flag, making the
public call reentrant while preserving `trace` and `ifault` semantics.

Farebrother is translated directly from the upstream Ruben/AS 204 code,
including its density output and `ifault` conventions.

Imhof's upstream C++ calls R's `Rdqagi` (QUADPACK) on `[0, infinity)`. The
Fortran port uses the supplied `r_mod::integrate` improper-integral facility,
as required by the reuse rule. `epsabs` and `epsrel` are mapped to the more
stringent of the two tolerances because `r_mod::integrate` exposes one
relative-tolerance control. The returned `abserr` is the helper's raw integral
error estimate, matching the scale returned by the upstream wrapper.

## Validation

The test suite includes numerical references produced directly from the
upstream C++ implementations for Davies and Farebrother. Imhof references were
independently evaluated from Imhof's published integral using high-accuracy
adaptive quadrature, and Liu references use the equivalent noncentral
chi-square formulation.

Manual compilation was performed with GNU Fortran 14.2 using Fortran 2018 and
implicit-interface errors enabled. FPM was not installed in the build
environment, so the manifest could not be exercised by the `fpm` executable.
The source/test/example layout follows normal FPM conventions.
