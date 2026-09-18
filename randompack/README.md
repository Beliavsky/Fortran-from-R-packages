# randompack - modern Fortran translation

This directory translates the computational core of the R package
`randompack` 0.1.10 to modern free-form Fortran with FPM. The upstream package
provides independent random-number-generator objects with multiple engines,
state management, continuous and discrete distributions, and multivariate
normal simulation.

The original source snapshot is retained under `upstream/`. See
`PROVENANCE.md`, `NOTICE.md`, `THIRD-PARTY-NOTICES`, and `LICENSE` before
redistributing the package.

## Build

Place this directory at the root of `Fortran-from-R-packages`, alongside the
shared linear-algebra package:

```text
Fortran-from-R-packages/
  randompack/
  rfortran-linalg/
```

Then run:

```text
cd randompack
fpm build
fpm test
fpm run --example basic_randompack
```

`rfortran-linalg` is used only by multivariate-normal simulation. No BLAS,
LAPACK, or translated dependency source is vendored here, and the manifest
contains no system `-lblas`/`-llapack` linkage.

## Main API

The public module is `randompack`.

```fortran
use randompack, only : dp, randompack_rng, randompack_rng_type

type(randompack_rng_type) :: rng
real(dp) :: x(100)
integer :: info

rng = randompack_rng('pcg64', seed=123)
call rng%normal(x, info=info)
```

The constructor defaults to `x256++simd`, matching the upstream R function.
The SIMD engine names are implemented by portable scalar emulation of the same
eight logical streams. They therefore preserve stream layout without requiring
AVX2, AVX-512, or NEON.

Available engine names, in the same table order as upstream, are:

```text
x256++simd
x256**simd
sfc64simd
x256++
x256**
x128+
xoro++
pcg64
sfc64
squares
philox
cwg128
ranlux++
chacha20
```

`randompack_engines(names, descriptions)` returns this metadata. It is
translated for convenience but is excluded from the computational-function
coverage denominator because it only reports engine metadata.

### Configuration and state

`randompack_rng_type` exposes type-bound procedures corresponding to the
upstream R6 methods:

- `seed`, `randomize`, `jump`, `advance`, and `set_state`;
- `pcg64_set_inc`, `cwg128_set_weyl`, `sfc64_set_abc`,
  `chacha_set_nonce`, `philox_set_key`, and `squares_set_key`;
- `duplicate`, `serialize`, `deserialize`, and `engine_name`.

Engine state in the Fortran API is represented as packed `integer(int64)`
words. Upstream R accepts vectors of 32-bit numeric words and packs adjacent
pairs into 64-bit state words. This representation difference is intentional;
the engine transitions themselves use the same bit patterns.

`serialize` writes a `randompack_snapshot` derived type rather than the
upstream opaque C raw-byte blob. A Fortran snapshot preserves the complete
Fortran stream position, buffered words, and configuration flags and can be
restored exactly by `deserialize`.

### Continuous distributions

The type-bound methods are:

```text
unif
normal
skew_normal
lognormal
gumbel
pareto
exp
gamma
chi2
beta
t
f
weibull
mvn
```

They retain the upstream distribution parameterizations and invalid-parameter
checks. `mvn` accepts a symmetric positive-semidefinite covariance matrix and
an optional mean vector. It uses `rfortran-linalg` for the ordinary Cholesky factorization and a
package-local translation of upstream `rp_dpstrf` for complete-pivoting
positive-semidefinite Cholesky fallback. Full-rank draws preserve the upstream
column-wise normal batching, while rank-deficient draws preserve its single
rank-sized normal block and matrix-product order.

### Discrete operations

The type-bound methods are:

```text
int
perm
sample
raw
```

Bounded integer generation follows the same Lemire-style rejection logic as
upstream, including 16-bit draws for small ranges and 32-bit draws for larger
R-style integer ranges. `perm` uses incremental Fisher-Yates. `sample` uses
Floyd sampling for small samples and reservoir sampling for larger samples.
The deterministic tests include an upstream-C-derived PCG64 seed-123 sequence
covering 64 bounded integers, a 25-element permutation, a 10-of-40 sample, and
the 31 raw bytes immediately following those operations. The Fortran
permutation/sample labels are one-based, so those two reference vectors are the
upstream zero-based values plus one.

## Reproducibility and compatibility

The deterministic seed expansion and the 14 engine transitions are translated
at the bit level. Tests contain fixed seed-123 reference vectors for all engines. The three SIMD engines use the upstream `2^253` xoshiro stream
separation and SFC64 `2^61` counter spacing. Uniform conversion preserves the
upstream 52-bit default and optional 53-bit mantissa conventions.

Normal and exponential sampling now use the upstream modified NumPy
256-strip Ziggurat tables, raw-word decoding, and vector consumption order.
For PCG64 seed 123, deterministic regression vectors for ten normal and ten
exponential draws match the upstream C implementation bit-for-bit. Ziggurat rejection/tail paths use a Fortran translation of the upstream
double-precision BSD/OpenLibm `log`, `log1p`, and `exp` routines. Forced-tail
PCG64 regression cases match upstream C bit-for-bit. Gamma now uses the upstream vector Marsaglia-Tsang proposal/rejection order,
including 53-bit acceptance uniforms, fixed OpenLibm acceptance logs, and the
shape-<1 vector boost convention. PCG64 seed-123 regression vectors for both
`shape >= 1` and `shape < 1` match upstream C bit-for-bit. Beta, Student-t,
and F now reproduce the upstream vector composition order as well; PCG64
seed-123 regression vectors of 17 draws for each method match upstream C
bit-for-bit. Lognormal, Gumbel, Pareto, Weibull, and skew-normal now also
follow the upstream transform and vector-consumption order. PCG64 seed-123
17-draw regression vectors for all five match upstream C bit-for-bit in
`bitexact=.true.` mode. Multivariate-normal generation now also follows the
upstream full-rank and singular execution paths. PCG64 seed-123 5-by-3
regressions for an SPD covariance and a rank-2 covariance requiring complete
pivoting match all 15 output values bit-for-bit, and the following six uniforms
match as well, pinning the post-MVN RNG position.

`randomize()` is also implementation-specific. Upstream first requests
cryptographic system entropy. Portable standard Fortran has no cryptographic
entropy API, so this translation combines `date_and_time`, `system_clock`, and
SplitMix64. Use `seed()` for reproducible work.

See `API-COVERAGE.md` for method-level details.

## Translation coverage

Package status: **substantial**.

**1 of 1 (100%)** exported computational R functions are mapped. This fraction
measures mapped computational functions under the stated coverage convention;
it does **not** mean complete R/R6 ABI compatibility, identical floating-point
streams for every distribution, or byte-compatible C serialization.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `randompack_rng` | `randompack_rng_mod` | `randompack_rng` plus `randompack_rng_type` type-bound computational methods | substantial |

The exported `randompack_engines()` function is an informational metadata
listing and is therefore excluded from the denominator. A translated
`randompack_engines` procedure is nevertheless provided.

The README coverage statement and `fpm.toml` `[extra.translation]` metadata use
the same 1-of-1 computational coverage basis.

## Validation

The deterministic test suite exercises engine vectors, stream jumps/advance,
state setters and restoration, exact upstream-derived normal/exponential,
gamma, beta, Student-t, F, lognormal, Gumbel, Pareto, Weibull, skew-normal,
and discrete vectors, distribution parameter/range behavior, discrete sampling,
and exact full-rank/pivoted-singular multivariate-normal regressions. `VALIDATION.md` records
which requested release checks were executable in the translation environment.

### Ziggurat compatibility

The maintained tests include upstream-C bit vectors for ordinary PCG64 normal/exponential draws and forced tail-path cases (normal seed 9071 and exponential seed 753). The tail/rejection path uses a Fortran translation of the double-precision BSD/OpenLibm `log`, `log1p`, and `exp` routines retained in the upstream source.
