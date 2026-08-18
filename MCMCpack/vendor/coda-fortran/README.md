# coda-fortran

`coda-fortran` is a modern Fortran/FPM translation of the computational core of the R package **coda 0.19-4.1**, "Output Analysis and Diagnostics for MCMC".

The original package is licensed under **GPL (>= 2)**. This translation is distributed under **GPL-2.0-or-later**. The original package metadata, authorship information, README, and changelog are retained under `original/`, and the GPL v2 text is in `COPYING`.

## Scope

Translated computational functionality:

- `mcmc` / `mcmc.list`-style chain containers, iteration metadata, pooling, and windowing
- autocorrelation matrices and diagonal autocorrelation summaries
- cross-correlation matrices
- rejection-rate diagnostic
- HPD intervals
- batch standard errors
- MCMC summaries: mean, SD, naive SE, time-series SE, and quantiles
- `spectrum0.ar`-style spectral density at zero using Yule-Walker AR fits with AIC order selection
- `spectrum0`-style low-frequency periodogram regression
- effective sample size
- Geweke diagnostic
- Gelman-Rubin diagnostic, upper confidence bound, optional transform, and multivariate PSRF
- Heidelberger-Welch diagnostic
- Cramer-von Mises CDF helper used by Heidelberger-Welch
- Raftery-Lewis diagnostic

Intentionally omitted because they are plotting/UI/I/O/R-object infrastructure rather than computational algorithms:

- trace, density, autocorrelation, cross-correlation, Gelman, Geweke, cumulative, and lattice plotting
- interactive `codamenu` infrastructure and graphics layout helpers
- JAGS/OpenBUGS/CODA file readers and format converters
- R S3 printing/coercion/data-frame plumbing that has no Fortran analogue

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example basic_diagnostics
```

The library has no external dependencies.

## Minimal use

```fortran
use coda

real(dp) :: draws(1000,2)
type(mcmc_chain) :: chain
type(geweke_result) :: gw
real(dp), allocatable :: ess(:)

! ... fill draws ...
chain = make_mcmc(draws)
ess = effective_size(chain)
gw = geweke_diag(chain)
```

## R-to-Fortran API mapping

| R coda | Fortran |
|---|---|
| `mcmc` | `make_mcmc` / `mcmc_chain` |
| `mcmc.list` | `make_mcmc_list` / `mcmc_list` |
| `window.mcmc` | `window_mcmc` |
| `autocorr` | `autocorr` |
| `autocorr.diag` | `autocorr_diag`, `autocorr_diag_list` |
| `crosscorr` | `crosscorr` |
| `rejectionRate` | `rejection_rate`, `rejection_rate_list` |
| `HPDinterval` | `hpd_interval` |
| `batchSE` | `batch_se`, `batch_se_list` |
| `spectrum0.ar` | `spectrum0_ar` |
| `spectrum0` | `spectrum0` |
| `effectiveSize` | `effective_size`, `effective_size_list` |
| `geweke.diag` | `geweke_diag` |
| `gelman.diag` | `gelman_diag` |
| `heidel.diag` | `heidel_diag` |
| `raftery.diag` | `raftery_diag` |
| `pcramer` | `cramer_cdf` |
| `summary.mcmc` | `summarize_mcmc` |
| `summary.mcmc.list` | `summarize_mcmc_list` |

## Numerical notes

The R package delegates several calculations to R's `stats` implementation. This translation supplies self-contained equivalents:

- AR order selection uses Yule-Walker/Levinson recursion and AIC.
- F quantiles use the regularized incomplete beta function and inversion.
- Normal quantiles use a high-accuracy rational approximation.
- The Cramer-von Mises helper evaluates the package's four-term formula and a self-contained approximation for `K_(1/4)`.
- `spectrum0` uses a direct DFT after the same optional batching to at most 200 observations; this keeps the implementation dependency-free and is appropriate for the small series used by that estimator.

The tests include fixed numerical checks for the Cramer-von Mises helper and Gelman-Rubin results, plus smoke/invariant tests for the remaining diagnostics.
