# RPESE modern Fortran

A modern Fortran 2018 translation of the computational core of the R package
`RPESE` (version 1.2.7), with an FPM package layout.

RPESE estimates standard errors for risk and performance measures using
influence functions, spectral long-run variance estimation, and bootstrap
methods. Plotting, `xts`/`zoo` containers, R object formatting, and R help-system
infrastructure are intentionally omitted.

## Implemented measures

- mean and robust mean
- standard deviation and semi-standard deviation
- value at risk and expected shortfall
- Sharpe, Sortino, and downside Sharpe ratios
- expected-shortfall and value-at-risk ratios
- Rachev ratio
- lower partial moments and Omega ratio

## Implemented standard-error methods

- `se_if_iid`: iid influence-function standard error
- `se_if_cor`: correlated influence-function standard error
- `se_if_cor_pw`: AR(1)-prewhitened correlated influence-function standard error
- `se_if_cor_adapt`: adaptive blend of the preceding two correlated methods
- `se_boot_iid`: iid nonparametric bootstrap
- `se_boot_cor`: fixed-length circular block bootstrap

The correlated influence-function methods fit a polynomial model to the
periodogram with the translated RPEGLMEN exponential or Gamma elastic-net
solver. RPEIF supplies all translated influence functions and robust cleaning.
Both dependencies are vendored for reproducible builds.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example rpese_demo
```

With GNU Make and GNU Fortran:

```text
make MODE=checked test examples
make clean
make MODE=optimized test examples
```

## Minimal example

```fortran
use rpese, only : dp, rpese_options, se_result, se_if_cor_pw, mean_se

real(dp) :: returns(100)
type(rpese_options) :: options
type(se_result) :: result

! Fill returns first.
options = rpese_options()
call mean_se(returns, result, se_if_cor_pw, options)
```

The generic `estimate_se` routine accepts a measure name and a method integer.
Named wrappers such as `es_se`, `var_se`, `sr_se`, and `omegaratio_se` mirror the
exported R entry points in Fortran-compatible spelling.

## Numerical notes

- The periodogram is computed by a self-contained discrete Fourier transform,
  avoiding an external FFT dependency. It is mathematically equivalent to the
  frequencies used by the R code but is O(n^2).
- AR(1) prewhitening uses least-squares estimation with an intercept rather than
  R's `arima` maximum-likelihood implementation.
- The block bootstrap uses circular fixed-length blocks, a deterministic and
  self-contained equivalent of the upstream fixed-block `tsboot` workflow.
- `source_compatibility=.true.` is the default. It preserves the upstream DSR
  point-estimate treatment of the risk-free rate and the upstream adaptive
  weighting rule for negative AR(1) coefficients. Setting it to `.false.` uses
  corrected definitions.

## License

The upstream RPESE package is GPL-2.0-or-later. This translation links the
GPL-3.0-or-later RPEIF translation, so the combined package is distributed under
GPL-3.0-or-later. The complete upstream source snapshot and license texts are
included.
