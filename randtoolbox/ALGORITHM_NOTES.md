# Algorithm notes

## Integer arithmetic

Unsigned 32-bit states are represented by nonnegative `integer(int64)` values
where convenient.  WELL uses `integer(int32)` bit patterns as in its standalone
port.  The LCG supports modulus 2^64 without signed-overflow assumptions by
forming the low 64 bits from 16-bit limbs.  For other large moduli, modular
multiplication uses overflow-safe repeated doubling.

## SFMT

The upstream implementation selects SSE2/Altivec/scalar code at compile time.
The Fortran port uses one scalar 128-bit-word recurrence represented as four
32-bit lanes.  This is portable and sequence-compatible with the scalar/native
upstream stream.  The 32 parameter sets for exponents 607 through 19937 and the
fixed sets for larger exponents are retained.

## Sobol

The 1111-dimensional initial direction-number table and primitive-polynomial
data are translated from `LowDiscrepancy-sobol-orig1111.c`.  The R package
currently warns that scrambling is disabled; the Fortran API therefore exposes
the effective pure sequence directly rather than a no-op scrambling argument.

## RNG state

The R package installs generators into R's global RNG callback machinery.
Modern Fortran instead uses explicit derived-type state, which makes multiple
independent streams possible and avoids process-global callbacks.

## Collision tests

The original `coll.test` accepts an arbitrary R function and repeatedly invokes
it.  The numerical parts are separated here: `collision_count` computes sample
collisions, and `collision_test_counts` evaluates the exact/Poisson reference
distribution and chi-square statistic from repeated collision counts.  Calling
an arbitrary user RNG is naturally done by the Fortran caller.
