# PeerPerformance Fortran

A modern Fortran 2018 and FPM translation of the computational routines in
`PeerPerformance` 2.4.0.

The upstream package implements the luck-corrected peer-performance framework
of Ardia and Boudt. This port works directly with dense real arrays and typed
results, has no external numerical dependencies, and preserves the upstream
GPL-2-or-later license and attribution.

## Features

- ordinary and factor-adjusted alpha estimates
- Sharpe and Cornish-Fisher modified-Sharpe ratios
- pairwise alpha, Sharpe, and modified-Sharpe tests
- classical and automatic Parzen-HAC covariance estimates
- IID and circular block bootstrap tests
- data-driven bootstrap block-length selection
- Storey-style null proportion estimation with finite-sample adjustment
- deterministic data-driven lambda selection
- positive, zero, and negative peer-performance ratios
- within-group and cross-group screening
- targeted fund screening
- rolling-window screening
- exposure-heterogeneity summaries
- IEEE-NaN complete-case handling

## Build

```sh
fpm build
fpm test
fpm run peerperformance_demo
fpm run --example basic_ratios
fpm run --example peer_screening
```

A direct GNU Fortran validation script is also supplied:

```sh
./scripts/build_validate.sh
```

On Windows with GNU Fortran:

```bat
scripts\build_validate.bat
```

## Basic use

```fortran
use peerperformance, only: dp, peer_control, screening_result, sharpe_screening

real(dp) :: returns(120,10)
type(peer_control) :: control
type(screening_result) :: result

control%has_lambda = .true.
control%lambda = 0.5_dp
control%min_obs = 30
call sharpe_screening(returns,control,result)

if (result%status /= 0) error stop trim(result%message)
```

For coefficient-level alpha screening:

```fortran
control%screen_beta = .true.
call alpha_screening(returns,control,result,factors=factors)
```

The first coefficient is the intercept/alpha. Subsequent coefficients correspond
to factor columns in their supplied order.

## Main result fields

`screening_result` contains:

- `estimate(coefficient,fund)`
- pairwise `difference`, `standard_error`, `tstat`, and `pvalue`
- `pipos`, `pizero`, and `pineg`
- selected `lambda`
- complete observation and peer counts

`test_result` contains the two estimates, their difference, standard error,
t-statistic, p-value, status, and diagnostic message.

## Numerical conventions

Missing observations are represented by IEEE NaN and removed pairwise. Standard
Sharpe ratios use the sample standard deviation with denominator `n-1`.
Modified Sharpe ratios use population central moments in the Cornish-Fisher
modified-VaR expression, matching the upstream implementation.

For HAC inference, this port provides an automatic Parzen-kernel long-run
covariance estimator. It is self-contained but is not intended to reproduce
bit-for-bit the exact `sandwich::vcovHAC` implementation used by R.

Bootstrap resampling follows the same IID/circular-block methodology as the R
package, but the random-number stream is Fortran-specific. A fixed `seed` gives
repeatable results within this implementation.

## Scope

The full numerical API exported by the upstream package is represented. R-only
S3 printing, plotting, data-frame conversion, parallel-cluster management,
formula/name handling, and loading of the bundled `hfdata.rdata` file are not
compiled. The original data and source are retained under `original/`.

See [COVERAGE.md](COVERAGE.md), [PORTING_NOTES.md](PORTING_NOTES.md), and
[VALIDATION.md](VALIDATION.md) for details.

## License

GPL-2.0-or-later. See `LICENSE` and `NOTICE`.

Upstream copyright:

Copyright (C) 2012-2023 David Ardia and Kris Boudt.
