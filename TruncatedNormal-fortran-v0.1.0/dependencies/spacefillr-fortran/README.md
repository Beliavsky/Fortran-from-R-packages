# spacefillr-fortran v0.1.0

Modern free-form Fortran/FPM port of the computational generators in
`spacefillr` 0.4.0.

Implemented generators:

- Faure-permuted Halton, up to 256 dimensions
- randomized-permutation Halton, up to 256 dimensions
- Joe-Kuo/PBRT Sobol, up to 1024 dimensions
- fast Owen-scrambled Sobol, up to 21201 dimensions
- progressive jittered (PJ)
- progressive multi-jittered (PMJ)
- blue-noise best-candidate PMJ
- PMJ(0,2)
- blue-noise best-candidate PMJ(0,2)

All public set routines return `real(real64)` matrices. Internally, PCG32,
Sobol, Owen scrambling, and the Halton generated kernels preserve the
upstream 32-bit and single-precision conversion semantics where those
semantics determine the sequence.

The library has no runtime dependencies. It uses `selected_int_kind(38)` for
exact modular 64-bit PCG arithmetic; gfortran 14 supports this integer kind.

## Random Halton portability

Upstream randomized Halton uses `std::shuffle`. The exact permutation produced
from a fixed engine state is not specified by the C++ standard and can differ
between standard-library implementations. This port instead specifies a
Fisher-Yates shuffle driven directly by the translated PCG32 generator. Thus
randomized Halton is reproducible across Fortran compilers implementing the
required integer kind, but it is intentionally not promised to be bit-for-bit
identical to every C++ `std::shuffle` implementation.

## FPM

The manifest uses free source form and disables both implicit typing and
implicit external interfaces. No BLAS/LAPACK or C/C++ library is required.
