# API reference

All real-valued calculations use `dp = kind(1.0d0)` from `frapo_kinds`.
Most routines that can fail accept an optional integer `status` argument.
Status constants are exported by module `frapo`.

## Result type

`type(portfolio_result)` contains:

- `weights(:)`: portfolio weights
- `drawdowns(:)`: realized drawdown path when applicable
- `objective`: optimized objective; terminal return for constrained-return
  drawdown portfolios and CDaR for `pmincdar`
- `terminal_return`: final cumulative portfolio return for drawdown portfolios
- `threshold`: empirical drawdown quantile for CDaR portfolios
- `risk_value`: variance objective, maximum drawdown, average drawdown, or CDaR,
  according to the portfolio type
- `status`, `iterations`, and `portfolio_type`

## Series routines

### `capser(y, minimum, maximum)`

Caps every value to `[minimum, maximum]`. Vector and matrix overloads are
provided. The descriptive name is `cap_series`.

### `returnseries(y, method, percentage, trim, compound, status)`

`method` is `returns_continuous` or `returns_discrete`. The default is
continuous returns. `percentage` defaults to true; `trim` and `compound`
default to false. Vector and matrix overloads are available.

The untrimmed first observation is IEEE NaN unless `compound` is true, in which
case it is zero.

### `returnconvert(y, direction, percentage)`

`direction` is `convert_cont_to_disc` or `convert_disc_to_cont`.

### Trend routines

- `trdbinary(y)`
- `trdbilson(y, exponent)`
- `trdes(y, lambda [, initial] [, status])`
- `trdhp(y, lambda [, status])`
- `trdsma(y, periods [, trim] [, status])`
- `trdwma(y, weights [, trim] [, status])`

The original names support vector and matrix data. Descriptive procedure names
begin with `trend_`.

## Risk and dependence routines

### `mrc(weights, covariance, percentage, status)`

Marginal contribution to portfolio volatility. Percentage output is the
default.

### `dr(weights, covariance, status)`

Diversification ratio.

### `cr(weights, covariance, status)`

Concentration ratio based on volatility-scaled positions.

### `rhow(weights, covariance, status)`

Volatility-weighted average pairwise correlation.

### `tdc(x, method, lower_tail, k, status)`

Nonparametric tail-dependence matrix. `method` is `tdc_empirical` or `tdc_evt`.
The default is empirical lower-tail dependence with
`k = floor(sqrt(number_of_observations))`. Average ranks are used for ties,
matching R's default `rank` behavior.

### `sqrm(matrix, status)`

Returns the real symmetric matrix square root computed by LAPACK eigendecomposition.

## Portfolio routines

Every routine returns `type(portfolio_result)`.

### `pgmv(returns [, percentage])`

Long-only fully invested global minimum-variance portfolio.

### `pmd(returns [, percentage])`

Most-diversified portfolio. It solves the correlation-space quadratic program,
rescales by individual volatilities, and renormalizes.

### `pmtd(returns [, method] [, k] [, percentage])`

Minimum-tail-dependence portfolio, followed by the original volatility
rescaling and normalization.

### `perc(covariance [, initial] [, percentage])`

Equal-risk-contribution portfolio using cyclic coordinate descent.

## Drawdown portfolio routines

Prices are passed as a `time x assets` matrix. All routines use cumulative
discrete returns from the first price row.

### `pmaxdd(price_data [, max_drawdown] [, soft_budget])`

Maximizes terminal cumulative return subject to a maximum-drawdown constraint.

### `pavedd(price_data [, average_drawdown] [, soft_budget])`

Maximizes terminal cumulative return subject to an average-drawdown constraint.

### `pcdar(price_data [, alpha] [, bound] [, soft_budget])`

Maximizes terminal cumulative return subject to a conditional-drawdown-at-risk
bound. Defaults are `alpha=0.95` and `bound=0.05`.

### `pmincdar(price_data [, alpha] [, soft_budget])`

Minimizes conditional drawdown at risk. The default is `alpha=0.95`.

When `soft_budget` is false, weights sum to one. When true, the sum is at most
one, matching the original package.
