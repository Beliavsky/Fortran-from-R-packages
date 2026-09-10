# dlm

Modern free-form Fortran translation of the computational core of the R package
**dlm 1.1-6.1**, with FPM packaging and deterministic tests.

The upstream package implements Bayesian and likelihood analysis for dynamic
linear models. This translation focuses on numerical model construction,
Kalman filtering/likelihood, smoothing, forecasting, simulation primitives, and
MCMC diagnostics, maximum-likelihood fitting, and adaptive rejection Metropolis
sampling. R plotting, printing, S3/list construction, time-series attributes,
and interactive progress displays are omitted; R closures and `...` forwarding
are represented by typed Fortran callback objects.

## Build

This directory is intended to live at the root of
`Beliavsky/Fortran-from-R-packages`, beside the shared packages it reuses:

```text
Fortran-from-R-packages/
  dlm/
  rfortran-core/
  rfortran-linalg/
  optimx/
```

Then run:

```text
cd dlm
fpm build
fpm test
fpm run --example local_level
fpm run --example arms_sampling
fpm run --example mle_variance
```

`rfortran-core` supplies the common `dp` kind. `rfortran-linalg` supplies the
checked dense linear algebra and uses the repository-standard pinned
`fortran-lapack` dependency, so `dlm` does not link system BLAS/LAPACK. The
existing top-level `optimx` translation supplies optimization for `dlm_mle`; no
optimizer implementation is copied into this package.

See `VALIDATION.md` for the exact validation environment, test results, source
audits, and the documented fact that FPM itself was unavailable in the build
sandbox.

## Main API

The public `dlm` module re-exports:

- `dlm_model`, `dlm_filter_result`, `dlm_smooth_result`, `dlm_forecast_result`
- `dlm_gibbs_dig_result`, `dlm_indicator`, `dlm_log_density`, `dlm_model_builder`
- `dlm_mod_reg`, `dlm_mod_poly`, `dlm_mod_seas`, `dlm_mod_trig`, `dlm_mod_arma`
- `ar_trans_pars`, `block_diag2`, `convex_bounds`, `dlm_add`, `dlm_sum`
- `dlm_filter`, `dlm_ll`, `dlm_smooth`, `dlm_forecast`, `dlm_residuals`
- `dlm_svd2var`, `dlm_bsample`, `dlm_gibbs_dig`, `rwishart`, `dlm_random_model`
- `dlm_mle`, `dlm_mle_result`, `arms`
- `mcmc_sd`, `mcmc_mean`, `erg_mean`, `seed_dlm_rng`

Observations are passed as a rank-2 array with time in rows. IEEE quiet NaNs
represent missing observation components. A time-varying model stores the same
R-style integer index maps (`JFF`, `JV`, `JGG`, `JW`): zero leaves a baseline
matrix entry unchanged and a positive value selects that column of `X` at the
current time row.

The filter and smoother expose full covariance arrays rather than R's lists of
SVD factors. This is a deliberate typed-Fortran API difference, not a claim of
bit-identical object compatibility.

## Translation coverage

Package status: **substantial**.

**25 of 25 (100%)** exported computational R functions have a meaningful
complete, substantial, or partial mapping. This fraction measures mapped
computational functions, **not** complete R interface compatibility or
bit-for-bit numerical identity.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bdiag` | `dlm` | `block_diag2` | partial |
| `convex.bounds` | `dlm` | `convex_bounds` | substantial |
| `arms` | `dlm` | `arms` | substantial |
| `ARtransPars` | `dlm` | `ar_trans_pars` | complete |
| `dlmModReg` | `dlm` | `dlm_mod_reg` | substantial |
| `dlmModPoly` | `dlm` | `dlm_mod_poly` | complete |
| `dlmModSeas` | `dlm` | `dlm_mod_seas` | complete |
| `dlmModTrig` | `dlm` | `dlm_mod_trig` | substantial |
| `dlmModARMA` | `dlm` | `dlm_mod_arma` | partial |
| `dlmSum` | `dlm` | `dlm_sum` | substantial |
| `%+%` | `dlm` | `dlm_sum` | substantial |
| `dlmMLE` | `dlm` | `dlm_mle` | substantial |
| `dlmLL` | `dlm` | `dlm_ll` | substantial |
| `dlmFilter` | `dlm` | `dlm_filter` | substantial |
| `dlmSmooth` | `dlm` | `dlm_smooth` | substantial |
| `dlmBSample` | `dlm` | `dlm_bsample` | substantial |
| `dlmGibbsDIG` | `dlm` | `dlm_gibbs_dig` | substantial |
| `rwishart` | `dlm` | `rwishart` | substantial |
| `dlmSvd2var` | `dlm` | `dlm_svd2var` | complete |
| `dlmForecast` | `dlm` | `dlm_forecast` | substantial |
| `dlmRandom` | `dlm` | `dlm_random_model` | partial |
| `mcmcSD` | `dlm` | `mcmc_sd` | substantial |
| `mcmcMean` | `dlm` | `mcmc_mean` | substantial |
| `mcmcMeans` | `dlm` | `mcmc_mean` | substantial |
| `ergMean` | `dlm` | `erg_mean` | complete |

All 25 functions in the exported-computational-function coverage basis now have
a mapping. `dlm_mle` uses a typed `dlm_model_builder` callback and the sibling
`optimx` package in place of R's arbitrary closure/`...`/`optim` list interface.
`arms` retains the upstream Gilks adaptive-rejection Metropolis envelope and
Metropolis correction, while typed density/support objects replace R closures.
Fortran intrinsic RNG sequences are not expected to match R bit-for-bit.

See `API_COVERAGE.md` and the `[extra.translation]` section of `fpm.toml` for
function-level notes. The README and manifest use the same 25/25 coverage basis.

## Licensing and provenance

The upstream package declares GPL version 2 or later. Maintained Fortran source
uses `SPDX-License-Identifier: GPL-2.0-or-later`; see `LICENSE`, `COPYING`,
`NOTICE`, and `PROVENANCE.md`. Review copies of the upstream source files used
for this translation are under `provenance/upstream/` and are not compiled.
