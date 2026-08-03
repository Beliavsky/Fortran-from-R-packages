# API

All real-valued inputs use `dp = kind(1.0d0)` from module `infoset`.

## Gross returns

```fortran
returns = g_ret(prices)
```

`prices` may be a rank-1 vector or a rank-2 matrix whose columns are assets.
Each output column is sorted in ascending order, matching the R function.

## Tail mixture

```fortran
call tail_mixture(y, shift, iteration, result [, control])
```

`result` is `tail_mixture_result` and contains the change point, ordered left
and right lognormal parameters, prior probability, type-I/type-II overlap
probabilities, log likelihood, convergence flag, and status.

## Information set

```fortran
call infoset_estimate(gross_returns, result [, control])
```

`result` is `information_set_result`. At most two valid change points below the
sample median are returned, matching the upstream `k < 3` loop.

## Overlapping windows

```fortran
windows = create_overlapping_windows(data, window_size, overlap)
```

The result has shape `(window_size, n_assets, n_windows)`, with consecutive
windows starting `overlap` observations apart. The upstream argument named
`ov` is therefore a step size, despite being described as an overlap.

## Left Risk

```fortran
call lr_cp(prices, window_size, overlap, result [, control])
```

`result%values` has shape `(n_assets, n_windows)`. The first change point is
estimated from the full price history and then used in each rolling window.
When no stable mixture split exists, the 10th percentile of full-sample gross
returns is used as a deterministic fallback and status is `infoset_no_split`.

## Portfolio construction

```fortran
call ptf_construction(prices, window_size, overlap, strategy, result, &
                      left_risk, penalty)
```

Supported strategies are:

- `M`: Markowitz covariance.
- `C_M`: Markowitz covariance with Left Risk linear term.
- `EDC`: extreme-downside-correlation matrix.
- `C_EDC`: EDC with Left Risk linear term.

`left_risk` is required for the combined strategies and has shape
`(n_assets, n_windows)`. The default penalty is `1.0e-4_dp`.

The QP imposes full investment, portfolio expected return equal to the
cross-sectional average expected return, and `0 <= weight <= 1`.

## Summary

```fortran
call summary_ptf(oos_returns, summary)
```

Returns count, minimum, first quartile, median, mean, third quartile, and
maximum in a `portfolio_summary` object.
