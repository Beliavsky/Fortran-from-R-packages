# API coverage and compatibility notes

## Scope

The coverage basis is **exported computational R functions**. It counts exported numerical functions and explicitly registered computational S3 methods. It excludes printing/formatting, character presentation, subsetting/replacement, dimensions/names, concatenation/binding/replication, matrix constructors, higher-order R callback interfaces (`apply`, `outer`), class predicates and interop helpers, metadata/version access, anonymous S4 glue, internal helpers, and `Math.bigz`/`Math.bigq` because those methods delegate their numerical work to the separate Rmpfr package.

Mapped coverage is **126 of 132 (95.5%)**; package status is **substantial**. The authoritative machine-readable mapping is in `fpm.toml`.

## Untranslated computational functions

- `is.na.bigz`, `is.na.bigq`: the current Fortran `bigz`/`bigq` types intentionally have no NA sentinel. Every representable value is a mathematical integer or normalized finite rational.
- `modulus`, `modulus.bigz`, `modulus<-`, `modulus<-.bigz`: the upstream R class can attach a modulus attribute to a `bigz` object and alter later arithmetic dispatch. The Fortran API instead makes modular arithmetic explicit with `bigz_modulo`, `bigz_powm`, and `bigz_invmod`; it does not attach mutable R-style object metadata.

## Material partial mappings

- Power/factorial/Fibonacci/combinatorial APIs use signed-64-bit loop/exponent/index arguments even though the values themselves are arbitrary precision.
- `isprime_z` uses deterministic Miller-Rabin for values fitting signed 64 bits and fixed bases for larger values; `factorize_z` uses Pollard rho. These are not GMP's exact algorithms/performance profile.
- `urand_bigz` uses a deterministic Park-Miller bit source. It is reproducible for a given seed in this translation, but not stream-compatible with GMP/R.
- `crossprod`/`tcrossprod` translate the self-product form; the optional second matrix is not part of these Fortran procedures.
- `solve.bigz`/`bigz_inverse` solve ordinary integer matrices exactly over the rationals. The R package's modulus-attached matrix solve mode is omitted.
- R recycling, vectorized mixed coercion, names/dimnames, S3/S4 dispatch, NA propagation, and R warning/error objects are not reproduced.

## Provenance

The upstream package metadata reports `gmp` 0.7-5.1, GPL (>= 2), by Antoine Lucas, Immanuel Scholz, Rainer Boehme, Sylvain Jasson, and Martin Maechler. Original sources and manuals are retained under `upstream/`; maintained Fortran files are a clean-room numerical translation guided by those sources and behavior.
