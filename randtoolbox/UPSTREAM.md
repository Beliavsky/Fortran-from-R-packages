# Upstream provenance

Translated from the supplied source tree for `randtoolbox` 2.0.5.

Major sources:

- `src/randtoolbox.c`, `src/congruRand.c`: package RNG and quasi-RNG kernels;
- `src/SFMT.c` and associated parameter headers: SFMT;
- `src/mt19937ar.c`: MT19937, 2002 initialization;
- `src/knuthTAOCP2002.c`: Knuth TAOCP generator;
- `src/LowDiscrepancy-sobol-orig1111.c`: Sobol direction data/recurrence;
- `src/LowDiscrepancy-halton.c`, `src/primes.h`: Halton/primes;
- `src/testrng.c` and `R/testRNG.R`: RNG test computations;
- `R/quasiRNG.R`, `R/pseudoRNG.R`, `R/binary-operator.R`, `R/sobol-basic.R`:
  public computational semantics.

The package depends on `rngWELL`; the previously translated native Fortran WELL
implementation is supplied by the canonical top-level `rngWELL` translation.
Its separate notices are retained in `LICENSES/`.
