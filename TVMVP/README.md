# TVMVP Modern Fortran

A self-contained modern Fortran translation of the computational routines in the R package **TVMVP 1.0.5** (repository directory `TV-MVP`). The library estimates a time-varying covariance matrix using kernel-weighted local principal components and adaptive POET residual covariance regularization, then applies the estimate to portfolio optimization and rolling evaluation.

## Included computations

- Silverman bandwidth and Epanechnikov/boundary/convolution kernels
- Local PCA at one time point or over a complete sample
- BIC-type factor-number selection
- Adaptive POET residual covariance thresholding
- Positive-definite covariance reconstruction
- Su-Wang constant-loading test with parametric bootstrap
- Mean, AR(1), MA(1), and ARMA(1,1) expected-return model selection
- Global minimum-variance, maximum-Sharpe, and return-constrained portfolios
- Expanding-window TV-MVP versus equal-weight evaluation
- Typed results, explicit status codes, and deterministic RNG seeding

Plotting, R6 classes, data frames, progress bars, `zoo`/tibble behavior, and formatted console output are omitted.

## Build

With FPM:

```text
fpm test
fpm run --example tvmvp_demo
```

With GNU Make:

```text
make MODE=checked test example
make MODE=optimized test example
```

The checked build uses Fortran 2018, warnings as errors, bounds/runtime checking, backtraces, and signaling-NaN initialization. The optimized build uses `-O3` and the same strict warning policy.

## Minimal use

```fortran
use tvmvp, only : dp, portfolio_prediction_result, predict_portfolio
real(dp) :: returns(100,20)
type(portfolio_prediction_result) :: result

! Fill returns, then:
call predict_portfolio(returns, 21, result, m=1, m0=5, &
                       compute_max_sharpe=.true., min_return=0.02_dp)
```

See `example/tvmvp_demo.f90` for a complete program.

## Source compatibility switches

Two upstream implementation details are preserved by default:

1. `boundary_kernel` normalizes boundary weights using the second index supplied by the R call site, not the local-PCA target index. Set `source_compatible_boundary=.false.` for target-index boundary normalization.
2. `expanding_tvmvp` refreshes its estimation sample and bandwidth only when the factor count is updated. Set `source_compatible_expanding=.false.` to refresh them at every rebalance.

The expected-return models use conditional sum-of-squares estimation rather than R's exact `stats::arima` likelihood. See `PORTING_NOTES.md`.

## License

MIT, matching the upstream package. The original source archive and license are retained under `upstream/` and `licenses/`.
