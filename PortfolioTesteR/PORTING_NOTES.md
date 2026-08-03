# Porting notes

## Array orientation

The Fortran port uses `time x asset` matrices. This keeps time slices contiguous
conceptually and makes the backtest interface explicit. Three-dimensional
feature and filter collections use `time x asset x feature/filter`. Parameter
combinations use `n_parameters x n_combinations`.

## Missing values

Missing numeric values are represented by IEEE quiet NaNs. Functions generally
ignore finite-data gaps where the upstream operation naturally permits it and
leave warmup or insufficient-data output as NaN. Logical selections are stored
as numeric zero/one matrices so they can be blended with partial regime weights.

## Dates and labels

The core routines do not carry dates, ticker names, sectors, or R attributes.
Callers retain those in parallel arrays. Group and sector membership is supplied
as integer IDs. `align_to_indices` performs exact integer-index matching with an
optional forward-fill mode.

## Indicators

The package's OHLC-dependent indicators are exposed in close-only forms because
the principal Fortran data representation is a price matrix. `calc_atr` is a
close-to-close true-range proxy. Rolling volatility supports standard deviation,
price range, robust MAD, mean absolute return, and downside deviation.

## Portfolio accounting

`run_backtest` is a share-and-cash engine. It sells before buying, can restrict
trades to whole shares, charges proportional basis-point costs, and records
executed rather than merely requested weights. Stop losses are evaluated at the
same observation frequency and price matrix supplied by the caller; no intraday
high/low path is inferred.

Weights are long-only. Negative and nonfinite requested weights are set to zero,
and rows above unit gross exposure are normalized. Unallocated weight remains
cash. The lightweight `portfolio_returns_from_weights` routine instead applies
lagged weights directly and is useful when detailed share accounting is not
needed.

## Risk-allocation methods

HRP uses native average-linkage clustering and recursive bisection. Equal-risk
contribution and maximum-diversification weights use self-contained iterative
methods rather than external R optimization packages. Results can differ by
small numerical amounts from a particular upstream solver or tie ordering.

## Machine learning

The self-contained model is pooled ordinary least squares or ridge regression.
Rolling predictions refit on each in-sample window. Random forests, gradient
boosting, and neural networks were not replaced with superficially similar
algorithms; callers may connect external implementations around the feature and
score APIs.

## Callback interfaces

`strategy_builder` is a Fortran procedure interface:

```fortran
subroutine builder(prices, params, weights, status)
  real(dp), intent(in) :: prices(:,:), params(:)
  real(dp), allocatable, intent(out) :: weights(:,:)
  integer, intent(out) :: status
end subroutine builder
```

This replaces R function objects and dynamic argument lists. Parameter-grid and
walk-forward routines call the procedure synchronously. Built-in strategies are
provided as examples and can be used directly.

## Numerical and language choices

- Real calculations use `dp = kind(1.0d0)`.
- Source is free-form Fortran 2018 with `implicit none` semantics enforced by
  the FPM manifest.
- Linear algebra is self-contained and intended for moderate educational and
  research problems, not as a replacement for optimized BLAS/LAPACK on very
  large panels.
- Random sample generation is deterministic when an explicit 64-bit seed is
  supplied.

## Licensing

The upstream package is MIT licensed. The translated source remains MIT
licensed, carries SPDX identifiers, and includes the complete upstream source
snapshot for provenance. No code from optional GPL or proprietary R dependencies
is bundled into the translated numerical modules.
