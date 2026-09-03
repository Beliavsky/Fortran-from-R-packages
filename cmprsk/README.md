# cmprsk

Modern free-form Fortran translation of the computational code in Robert Gray's R package **cmprsk** 2.2-12.

The translation focuses on competing-risks estimation and regression. R plotting, printing, S3 dispatch, formula/data-frame handling, and NA-row orchestration are intentionally omitted.

## Implemented numerical API

The public module is `cmprsk` and re-exports `dp` from `rfortran-core`.

- `fit_cuminc` — cumulative-incidence curves for every noncensoring cause and group, plus stratified Gray tests.
- `cumulative_incidence` — direct translation of the upstream `cinc` kernel, including its variance estimator.
- `gray_test` — direct translation of `crstm`/`crst` for stratified K-sample tests with arbitrary `rho`.
- `cuminc_timepoints`, `curve_timepoints`, `timepoint_indices` — numerical counterpart of `timepoints`/`tpoi`.
- `fit_crr` — Fine-Gray subdistribution-hazard regression with fixed effects, time-varying effects, censoring groups, Newton fitting, Armijo step reduction, robust sandwich covariance, score residuals, and baseline hazard jumps.
- `predict_crr` — cumulative-incidence predictions at supplied covariate profiles.
- `summarize_crr` — coefficient standard errors, z tests, normal p-values, relative risks, confidence intervals, and the pseudo-likelihood-ratio statistic.

`fit_crr` accepts the time-function matrix directly. This replaces the R callback `tf` with an explicit numerical interface: its rows correspond to ascending distinct failures of the modeled cause.

## Dependencies

This directory is intended to live beside the repository's shared packages:

```toml
rfortran-core   = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

`rfortran-core` supplies `dp` and R-compatible distribution functions. `rfortran-linalg` supplies checked dense solves/inverses and uses the repository's pinned pure-Fortran LAPACK dependency. No BLAS, LAPACK, or translated dependency source is copied into this package.

## Build and test

From this directory in the repository root:

```text
fpm build
fpm test
fpm run --example competing_risks_example
```

The maintained source does not require `-ffast-math` and should not be compiled with finite-math assumptions, because input NaNs are explicitly rejected with `ieee_is_nan`.

## Input conventions

Times may be unsorted in the high-level `fit_cuminc` and `fit_crr` interfaces; they are stably sorted internally. Integer cause/group/stratum labels may be arbitrary. Direct low-level kernels document when sorted input or consecutive group codes are required.

Unlike the R wrappers, the Fortran interface does not silently omit rows containing missing values. NaN numerical inputs are rejected, so callers can make their desired missing-row policy explicit before fitting.

See `API_COVERAGE.md`, `PROVENANCE.md`, and `VERIFICATION.md` for detailed parity and verification information.
