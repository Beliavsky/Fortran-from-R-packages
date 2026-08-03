# Porting notes

## Quantiles and ES tails

Plain historical simulation uses the R default type-7 sample quantile. As in
the upstream package, ES averages losses **strictly greater** than VaR. A sample
with no strict exceedance returns IEEE NaN and status `quarks_empty_tail`.

Age weighting preserves the upstream ordering, normalized exponential weights,
cumulative-probability interpolation, and weighted strict-tail mean.

## Volatility models

EWMA is a direct translation. The initial variance is the ordinary sample
variance of the complete input series, followed by

`h(t) = lambda*h(t-1) + (1-lambda)*x(t-1)^2`.

The GARCH route uses the included GPL-3.0-only `rugarch` Fortran port and fits a
Gaussian zero-mean sGARCH(1,1). If fitting fails, the routine returns an EWMA
forecast, sets `used_fallback=.true.`, and reports `quarks_fit_failed`.

## Filtered simulation RNG

R uses `sample()`. The Fortran API optionally accepts a stateful xorshift RNG,
allowing exact reproducibility without changing Fortran's global random state.
When no state is supplied, `random_number` is used. The Fortran interface uses a practical default bootstrap size of 10000; the R function defaults to `NULL` and requires the caller to supply it.

## Rolling scale smoothing

The original package stores an internal `smooth.help` object and delegates to
the R package `smoots`. This object is not ordinary exported source. The
Fortran `smooth_lpr` and `smooth_auto` modes use Gaussian-kernel smoothing of
squared returns and a one-step endpoint extrapolation. These modes are useful
adaptations, not numerical reproductions of `smoots`. `smooth_none` is exact.

## Coverage-test compatibility

The upstream conditional-likelihood branch uses `n10` as the exponent of
`p01` when `n11 > 0`; the standard Christoffersen expression uses `n01`.
`cvgtest` defaults to the upstream expression. Pass
`upstream_formula=.false.` to use the standard corrected transition count.

## Portfolio exact-return compatibility

For a time-varying weight matrix and `approxim=0`, upstream R subtracts one from
each asset-level product before the row sum. `plop_time_varying` preserves that
behavior by default. Pass `upstream_exact_mode=.false.` for the conventional
weighted simple return `sum((exp(x)-1)*w)`.

## Omitted infrastructure

Plotting, S3 printing, Shiny, Yahoo Finance downloading, progress bars,
`xts`/date classes, and packaged `.rda` datasets are not part of the Fortran API.
The complete upstream package is retained under `reference/` for provenance.
