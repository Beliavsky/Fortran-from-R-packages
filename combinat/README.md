# combinat

Modern free-form Fortran translation of the computational routines in the R package `combinat` 0.0-8.

The upstream package was written primarily by Scott Chasalow, with multinomial-simulation code derived from contributions by John Wallace and Alan Zaslavsky. The original R sources and documentation are retained under `upstream/` for provenance. This translation is unofficial and is not endorsed by the upstream authors, CRAN, or the R Foundation.

## Build

From a top-level `combinat` directory in the repository:

```text
fpm build
fpm test
```

No BLAS, LAPACK, ARPACK, C library, or translated R-package dependency is required. All maintained Fortran source is free form and uses the single `dp` kind exported by `combinat_kinds` and re-exported through `combinat`.

Examples can be run with:

```text
fpm run --example basic_usage
fpm run --example multinomial_usage
```

## Public API

Import the facade module:

```fortran
use combinat
```

The principal public operations are:

- `combn`, `combn2`: combinations and unordered pairs for indices, integer values, or real values.
- `permn`: all permutations in the upstream minimal-change order.
- `xsimplex`, `nsimplex`: simplex-lattice generation and counts.
- `hcube`: hypercuboid lattice generation with optional recycled scaling and translation.
- `ncm`, `ncm_vec`, `fact`, `logfact`: combinatorial and gamma-based helpers.
- `dmnom`: multinomial probability mass.
- `rmultinomial`, `rmultz2`: multinomial simulation using an explicit `combinat_rng_state`.
- `x2u`: expansion of bin counts to bin indices or integer/real labels.

R callbacks passed through `fun`, R list/array simplification, S-language missing-argument dispatch, and arbitrary R object types are interface behavior rather than numerical kernels and are not reproduced. Random simulation uses a deterministic Park-Miller generator exposed through `rng_seed`; it is distributionally equivalent in purpose but does not reproduce R RNG streams.

## Translation coverage

Package status: **substantial**.

**13 of 13 (100%)** exported computational R functions have meaningful Fortran mappings. This fraction measures mapped computational functions, not complete compatibility with every R calling convention, object type, callback, warning, or RNG behavior.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `combn` | `combinat_api` | `combn_indices`, `combn_integer`, `combn_real` | substantial |
| `combn2` | `combinat_api` | `combn2_indices`, `combn2_integer`, `combn2_real` | substantial |
| `dmnom` | `combinat_api` | `dmnom` | complete |
| `fact` | `combinat_api` | `fact` | complete |
| `hcube` | `combinat_api` | `hcube` | complete |
| `logfact` | `combinat_api` | `logfact` | complete |
| `nCm` | `combinat_api` | `ncm`, `ncm_vec` | substantial |
| `nsimplex` | `combinat_api` | `nsimplex`, `nsimplex_vec` | complete |
| `permn` | `combinat_api` | `permn_indices`, `permn_integer`, `permn_real` | substantial |
| `rmultinomial` | `combinat_api` | `rmultinomial` | substantial |
| `rmultz2` | `combinat_api` | `rmultz2` | substantial |
| `x2u` | `combinat_api` | `x2u_indices`, `x2u_integer`, `x2u_real` | substantial |
| `xsimplex` | `combinat_api` | `xsimplex` | substantial |

See `API_COVERAGE.md` for the interface differences and `NOTICE.md` for upstream provenance and licensing.

## Numerical and API notes

`ncm` follows the original `nCm` branching strategy for integer-valued `m`, including negative-`n` reflection and the finite-product case for positive noninteger `n < m`. `ncm_vec` and `nsimplex_vec` reproduce R-style recycling for nonempty vectors.

Allocation-producing routines return a status code when requested. `combinat_success` is zero; `combinat_invalid_argument` and `combinat_size_overflow` report invalid dimensions or impractical output sizes. This replaces R's `stop()`/warning behavior with Fortran-friendly explicit status handling.

The source does not rely on `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or similar assumptions. Explicit NaN tests use the intrinsic IEEE arithmetic facilities.
