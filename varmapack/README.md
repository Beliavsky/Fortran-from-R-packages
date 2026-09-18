# varmapack

Modern free-form Fortran translation of the computational core of the R package
`varmapack` 0.1.1 by Kristján Jónasson.

The package provides Gaussian VAR, VMA, VARMA, and VARMAX model construction;
burn-in-free stationary simulation; simulation conditional on startup values;
sample and theoretical autocovariances; PSI and orthogonalized impulse
responses; AR/MA spectral radii; and the upstream testcase collection.

## Dependencies

Place this directory at the same repository level as:

- `randompack`
- `rfortran-linalg`

The FPM manifest uses sibling path dependencies and does not vendor their
sources. No system BLAS or LAPACK link flags are used by this package.

```toml
[dependencies]
randompack = { path = "../randompack" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

## Build

From this directory in the target repository:

```text
fpm build
fpm test
fpm run --example basic_varmapack
```

## Basic use

```fortran
use varmapack, only : dp, varmapack_model, varmapack_model_type
use randompack, only : randompack_rng, randompack_rng_type

type(varmapack_model_type) :: model
type(randompack_rng_type) :: rng
real(dp), allocatable :: a(:, :, :), x(:, :, :), e(:, :, :)
integer :: info

allocate(a(1, 1, 1))
a = 0.5_dp
model = varmapack_model(reshape([1.0_dp], [1, 1]), a=a, info=info)
rng = randompack_rng(seed=123)
call model%sim(100, x, e, info, nrep=3, rng=rng)
```

Coefficient arrays use shape `(r,r,p)` for AR terms, `(r,r,q)` for MA terms,
and `(r,d,s)` for exogenous terms. Simulated series and shocks have shape
`(r,n,nrep)`.

## Public numerical API

The public facade module `varmapack` exports:

- `dp`
- `varmapack_model` and `varmapack_model_type`
- `varmapack_autocov`
- `varmapack_cov2corr`
- `varmapack_testcase`
- `varmapack_testcases`

A model exposes type-bound methods `sim`, `acvf`, `psi`, `irf`, `specrad`, and
`ma_specrad`.

`varmapack_autocov` uses ML normalization by default. Pass
`corrected=.true.` for the upstream `norm="C"` convention.

## Simulation

For a stationary VARMA model without supplied startup observations, simulation
begins from the exact stationary Gaussian state implied by the model rather
than using a burn-in period. When startup observations are supplied, the
startup shocks are drawn from their model-implied conditional Gaussian
distribution. Nonstationary VARMA models require supplied startup observations.

VARMAX simulation accepts exogenous input arrays with shape `(d,n,1)` for a
shared path or `(d,n,nrep)` for replicate-specific paths. Startup arrays use the
same one-path-or-`nrep` convention.

## Translation coverage

**Package status: substantial.**

**5 of 5 (100%)** exported computational R functions have meaningful Fortran
mappings. This percentage measures mapped computational functions; it does
**not** mean complete compatibility with R interfaces, R6 behavior, internal
linear-algebra paths, or floating-point/RNG ordering.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `varmapack_autocov` | `varmapack_analysis` | `varmapack_autocov` | complete |
| `varmapack_cov2corr` | `varmapack_analysis` | `varmapack_cov2corr` | complete |
| `varmapack_model` | `varmapack_model_mod` | `make_varmapack_model`, `model_sim`, `model_acvf`, `model_psi`, `model_irf`, `model_specrad`, `model_ma_specrad` | substantial |
| `varmapack_testcase` | `varmapack_testcases_mod` | `testcase_by_name`, `testcase_by_index` | substantial |
| `varmapack_testcases` | `varmapack_testcases_mod` | `varmapack_testcases` | substantial |

The R6 model's computational methods are included in the `varmapack_model`
mapping rather than counted as separate exported R functions.

### Main compatibility differences

- The API is typed Fortran rather than R6. Status failures use integer `info`
  codes instead of R conditions.
- `model%sim` always returns both the series and innovations as output
  arguments instead of varying the R return object with `return_shocks`.
- The stationary covariance setup uses an augmented-state discrete Lyapunov
  fixed-point iteration rather than upstream's VYW/SLICOT selection logic.
- Singular innovation covariances are supported; an eigensquare-root fallback
  is used for singular IRFs if ordinary Cholesky is unavailable.
- The sibling Fortran `randompack` supplies RNGs. Exact R/C RNG consumption is
  not claimed for every startup-factorization path.
- Printing, R data frames, R6 reflection, and other presentation behavior are
  not translated.

See `API-COVERAGE.md`, `PROVENANCE.md`, and `VALIDATION.md` for details.

## License and provenance

The upstream package is MIT licensed, copyright 2026 Kristján Jónasson. The
MIT license, original upstream metadata, selected source references, tests, and
third-party notice are retained. See `LICENSE`, `NOTICE.md`, and
`PROVENANCE.md`.
