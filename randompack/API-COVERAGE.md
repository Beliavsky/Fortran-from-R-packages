# API coverage

## Coverage basis

`NAMESPACE` exports `randompack_rng` and `randompack_engines`.
`randompack_rng` constructs the R6 object through which all computational
methods are exposed, so it is the single exported computational R function in
the formal coverage denominator. `randompack_engines` only reports static
engine names/descriptions and is excluded as informational metadata.

Formal function coverage is therefore **1 of 1 (100%)**, with package status
**substantial**. The method-level table below is more informative than that
single exported-function count.

## R6 method mapping

| Upstream method | Fortran type-bound procedure | Status | Main difference |
| --- | --- | --- | --- |
| `seed` | `seed` | substantial | Packed Fortran state representation; seed mixing/engine state is bit-level translated. |
| `randomize` | `randomize` | substantial | Uses portable clock/date entropy rather than OS cryptographic entropy. |
| `jump` | `jump` | substantial | Supported jump powers and engine transitions translated. |
| `advance` | `advance` | substantial | PCG64 128-bit delta is supplied as two packed 64-bit words. |
| `set_state` | `set_state` | substantial | Fortran accepts packed 64-bit words rather than R 32-bit numeric words. |
| `pcg64_set_inc` | `pcg64_set_inc` | substantial | Packed 128-bit increment. |
| `cwg128_set_weyl` | `cwg128_set_weyl` | substantial | Includes the upstream recommended 96-state warm-up. |
| `sfc64_set_abc` | `sfc64_set_abc` | substantial | Includes the upstream 18-state skip. |
| `chacha_set_nonce` | `chacha_set_nonce` | substantial | Three 32-bit nonce words are stored in `int64` arguments. |
| `philox_set_key` | `philox_set_key` | substantial | Packed 64-bit key words. |
| `squares_set_key` | `squares_set_key` | substantial | 64-bit key word. |
| `duplicate` | `duplicate` | substantial | Copies complete Fortran object/cache state. |
| `serialize` | `serialize` | substantial | Returns `randompack_snapshot`, not the upstream opaque byte blob. |
| `deserialize` | `deserialize` | substantial | Restores Fortran snapshots only. |
| `unif` | `unif` | substantial | 52/53-bit conversion and engine consumption translated. |
| `normal` | `normal` | substantial | Upstream Ziggurat tables/word ordering plus BSD/OpenLibm rejection/tail math translated; ordinary and forced-tail PCG64 reference vectors match upstream C bit-for-bit. |
| `skew_normal` | `skew_normal` | substantial | Upstream whole-vector first-normal draw followed by 128-value second-normal batches is translated; tested PCG64 vector matches bit-for-bit. |
| `lognormal` | `lognormal` | substantial | Upstream normal-vector, shift/scale, and BSD/OpenLibm exponentiation order is translated; tested PCG64 vector matches bit-for-bit. |
| `gumbel` | `gumbel` | substantial | Upstream vector-uniform and two-stage BSD/OpenLibm logarithm transform order is translated; tested PCG64 vector matches bit-for-bit. |
| `pareto` | `pareto` | substantial | Built from the upstream Ziggurat exponential stream followed by divide/exp/scale ordering; tested PCG64 vector matches bit-for-bit. |
| `exp` | `exp` | substantial | Upstream Ziggurat tables/word ordering plus BSD/OpenLibm rejection/tail math translated; ordinary and forced-tail PCG64 reference vectors match upstream C bit-for-bit. |
| `gamma` | `gamma` | substantial | Upstream vector Marsaglia-Tsang ordering, 53-bit acceptance uniforms, OpenLibm acceptance logs, and shape-<1 boost ordering are translated; tested PCG64 vectors match bit-for-bit. |
| `chi2` | `chi2` | substantial | Uses the translated upstream-order gamma construction. |
| `beta` | `beta` | substantial | Upstream whole-vector first-gamma then 128-value second-gamma batching and open-interval clamp translated; tested PCG64 vector matches bit-for-bit. |
| `t` | `t` | substantial | Upstream whole-vector normal draw followed by 128-value gamma batches translated; tested PCG64 vector matches bit-for-bit. |
| `f` | `f` | substantial | Upstream whole-vector numerator gamma followed by 128-value denominator gamma batches translated; tested PCG64 vector matches bit-for-bit. |
| `weibull` | `weibull` | substantial | Built from the upstream Ziggurat exponential stream with the same shape-1/shape-2 special cases and log/scale/exp ordering; tested PCG64 vector matches bit-for-bit. |
| `mvn` | `mvn` | substantial | Ordinary Cholesky uses shared `rfortran-linalg`; fallback translates upstream complete-pivoting PSD Cholesky and preserves full-rank/singular normal batching. Tested SPD and rank-2 PCG64 cases, including post-call stream position, match upstream C bit-for-bit. |
| `int` | `int` | substantial | Same unbiased small/large integer rejection strategy. |
| `perm` | `perm` | substantial | Same incremental Fisher-Yates construction. |
| `sample` | `sample` | substantial | Same Floyd/reservoir selection threshold; Floyd membership uses a linear search rather than a hash table. |
| `raw` | `raw` | substantial | Same buffered engine bytes on little-endian bit ordering. |

## Engine mapping

All upstream engine identifiers are accepted:

| Engine | Translation |
| --- | --- |
| `x256++simd` | Eight xoshiro256++ streams, scalar-emulated with upstream `2^253` separation. |
| `x256**simd` | Eight xoshiro256** streams, scalar-emulated with upstream `2^253` separation. |
| `sfc64simd` | Eight SFC64 streams with upstream `2^61` counter offsets. |
| `x256++` | xoshiro256++ bit-level transition. |
| `x256**` | xoshiro256** bit-level transition. |
| `x128+` | xorshift128+ bit-level transition. |
| `xoro++` | xoroshiro128++ bit-level transition. |
| `pcg64` | PCG64-DXSM with arbitrary advance and power-of-two jump. |
| `sfc64` | SFC64 bit-level transition. |
| `squares` | Squares64 bit-level transition. |
| `philox` | Philox-4x64-10 bit-level transition. |
| `cwg128` | CWG128 using portable 128-bit emulation. |
| `ranlux++` | ranlux++ 576-bit modular transition/jump using 64-bit limbs. |
| `chacha20` | Portable ChaCha20 block function and counter/nonce layout. |

The fixed-vector tests verify deterministic seed expansion and the first four
64-bit outputs of every engine. The SIMD tests additionally pin their lane
initialization and deterministic output sequences. PCG64 seed-123 tests also
pin ten upstream C normal draws and ten upstream C exponential draws bit-for-bit.
