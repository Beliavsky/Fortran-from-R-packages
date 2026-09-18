# Provenance

## Upstream

- R package: `randompack`
- Version translated: 0.1.10
- Upstream author and maintainer: Kristjan Jonasson, University of Iceland
- Upstream project: `https://github.com/jonasson2/randompack`
- CRAN publication recorded by the supplied source: 2026-08-24
- Upstream snapshot: retained verbatim under `upstream/`, except that no
  generated build products are included.

The original package is primarily an R6 interface over C implementations of
RNG engines, buffering/state management, probability distributions, and
multivariate-normal simulation. The translation follows the R-visible
computational behavior and the C implementation where that behavior is
defined.

## Translated numerical material

The maintained Fortran sources implement:

- `seed_seq_fe128`-style deterministic seed expansion;
- xoshiro256++, xoshiro256**, xorshift128+, xoroshiro128++, PCG64-DXSM,
  SFC64, Squares64, Philox-4x64, CWG128, ChaCha20, and ranlux++;
- portable scalar emulation of the upstream xoshiro256++ SIMD,
  xoshiro256** SIMD, and SFC64 SIMD stream layouts;
- power-of-two jumps, PCG64 arbitrary advance, and engine-specific state/key
  setters;
- 52-bit and 53-bit uniform conversion, unbiased bounded integer generation,
  permutations, sampling without replacement, and raw-byte draws;
- normal, skew-normal, lognormal, Gumbel, Pareto, exponential, gamma,
  chi-square, beta, Student-t, F, Weibull, and multivariate-normal draws.

## Deliberate implementation differences

The upstream modified NumPy Ziggurat normal/exponential tables, proposal
decoding, gap-bound tests, and vector prefetch/reverse-processing order are now
translated. The retained constants originate in the upstream source snapshot
and therefore remain under their applicable upstream third-party notices.
PCG64 seed-123 normal and exponential regression vectors match upstream C
bit-for-bit. The double-precision BSD/OpenLibm `log`, `log1p`, and `exp` routines used by
Ziggurat rejection/tail paths are translated in `src/randompack_openlibm.f90`.
Forced-tail PCG64 cases match upstream C bit-for-bit. Gamma now follows the
upstream vector Marsaglia-Tsang proposal/rejection order, fixed OpenLibm
acceptance logs, 53-bit acceptance uniforms, and shape-<1 boost ordering.
PCG64 reference vectors for both shape regimes match upstream C bit-for-bit.
Beta, Student-t, and F now also follow the upstream vector batching/order;
17-value PCG64 seed-123 references for each match upstream C bit-for-bit.
Lognormal, Gumbel, Pareto, Weibull, and skew-normal now reproduce the upstream
transform and vector-consumption order as well; 17-value PCG64 seed-123
reference vectors for all five match upstream C bit-for-bit in bitexact mode.

The upstream serialized state is an opaque C byte blob. The Fortran API uses a
strongly typed `randompack_snapshot` containing the complete Fortran engine,
buffer, and distribution-cache state. It restores Fortran streams exactly but
is not byte-compatible with the upstream C serialization format.

Upstream `randomize()` first requests cryptographic system entropy and has a
fallback. Standard Fortran has no portable cryptographic entropy API, so the
translation derives non-reproducible initialization from `date_and_time`,
`system_clock`, and SplitMix64. Deterministic `seed()` behavior is the route for
reproducible calculations.

For multivariate normal simulation, the upstream code uses ordinary LAPACK
Cholesky and then `rp_dpstrf`, its retained complete-pivoting PSD Cholesky
fallback. The translation calls shared `rfortran-linalg` for ordinary Cholesky
and translates the lower-triangle `rp_dpstf2` fallback directly in
`src/randompack_rng_mod.f90`, using the upstream fixed tolerance `1.0e-14`. It
then undoes row pivoting exactly as upstream before drawing. PCG64 seed-123
full-rank and rank-2 singular 5-by-3 regression matrices, plus the six uniforms
immediately following each MVN call, match the upstream C implementation
bit-for-bit in the tested cases.
