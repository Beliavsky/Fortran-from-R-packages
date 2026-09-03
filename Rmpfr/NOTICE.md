# NOTICE and provenance

## Upstream package

This work translates computational code from the R package **Rmpfr**, upstream
snapshot version **1.1-2**, dated 2025-10-21 in the supplied `DESCRIPTION`.
The supplied source identifies the following contributors:

- Martin Maechler — author and maintainer
- Richard M. Heiberger — contributor to hexadecimal/binary/decimal formatting
- John C. Nash — contributor of `hjkMpfr()` and origin of `unirootR()`
- Hans W. Borchers — contributor of `optimizeR(..., "GoldenRatio")` and origin
  of the Hooke-Jeeves code used by Rmpfr
- Mikael Jagan — contributor to safer conversion/configuration code

The upstream package is licensed **GPL (>= 2)**.  The Fortran translation is
therefore distributed as **GPL-2.0-or-later**; `LICENSE` contains GPL version 2.
The original `DESCRIPTION`, `README.md`, `NEWS.Rd`, and `ChangeLog` are retained
under `upstream/` for provenance.

## Native libraries

Rmpfr's numerical representation and arithmetic are defined by GNU MPFR and
GMP.  This translation uses those libraries through `ISO_C_BINDING` and does
not copy or vendor their sources.  Their own licenses remain applicable to the
installed libraries.  Upstream Rmpfr lists MPFR >= 3.2.0 and GMP >= 4.2.3 as
system requirements.

## Algorithm provenance retained from upstream comments

- `integrateR`: Bauer (1961), Algorithm 60, Romberg integration.
- Brent optimizer: based on the `localmin` algorithm described in Richard
  Brent, *Algorithms for Minimization without Derivatives* (1973), as recorded
  in upstream `optimizers.R`.
- Golden-ratio optimizer: contributed through Hans W. Borchers/pracma lineage,
  as recorded by upstream Rmpfr.
- `unirootR`: John C. Nash's R translation lineage of R's `zeroin` routine,
  as documented in upstream `unirootR.R`.
- Hooke-Jeeves optimization: John C. Nash/Hans W. Borchers contribution history
  as documented in upstream `hjk.R`.

## Dependency review for Fortran-from-R-packages

Before implementation, the top-level `Fortran-from-R-packages` repository was
checked for existing Rmpfr, MPFR, or GMP translations/shared modules.  None was
listed.  Consequently this package does not copy a sibling implementation and
instead binds to the same native MPFR/GMP libraries required by upstream.
No BLAS, LAPACK, ARPACK, `rfortran-compat`, or translated R-package dependency
source is included.

All ordinary Fortran real-valued code uses the single package kind `dp = real64`.
The low-level `ISO_C_BINDING` interface uses `real(c_double)` only for the three
MPFR C ABI entry points whose signatures literally contain a C `double`
(`mpfr_set_d`, `mpfr_get_d`, and `mpfr_cmp_d`); this is an interoperability
boundary, not a second numerical working precision.

## Intentional translation differences and corrections

1. **`tanpi` half-integers.**  The supplied upstream `src/Ops.c` returns `+1`
   or `-1` at reduced arguments `+/-0.5`.  Mathematically, `tan(pi*x)` has a
   pole there.  The Fortran translation deliberately returns signed infinity
   and tests this behavior.  This is an explicit numerical correction, not a
   claim of bit-for-bit parity with that upstream branch.
2. **Hooke-Jeeves direction order.**  Upstream `.hjexplore()` traverses the
   coordinate directions in a fresh `sample.int()` order.  The Fortran version
   scans coordinates in deterministic increasing index order.  The search
   steps and pattern-move logic are otherwise retained; deterministic ordering
   makes tests reproducible without importing R's RNG.
3. **qnorm bracketing.**  Upstream `qnormI()` uses specialized asymptotic
   starting intervals.  The Fortran version uses precision-preserving expanding
   brackets followed by the translated Brent root logic.  It favors a smaller,
   dependency-free core over reproducing R's double-precision `qnorm` helper.
4. **Array ownership API.**  Matrix and cumulative-array computations are
   output subroutines rather than allocatable-array function returns.  This
   avoids shallow-copy hazards for finalizable values containing MPFR storage
   across Fortran compilers while preserving the numerical operation.
5. R-specific S4/S3 dispatch, R vector recycling, printing/formatting methods,
   R serialization, and R `gmp` object conversions are interface functionality
   and are intentionally omitted.
