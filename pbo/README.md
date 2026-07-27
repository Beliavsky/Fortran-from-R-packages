# pbo-fortran

A modern Fortran/FPM translation of the computational portions of the R package
`pbo` 1.3.5, Probability of Backtest Overfitting.

The library implements combinatorially symmetric cross-validation (CSCV), rank
logits, probability of backtest overfitting, selected-strategy degradation and
loss diagnostics, selection frequencies, and stochastic-dominance curve data.
It is self-contained and has no external numerical dependencies.

## Build

```text
fpm build
fpm test
fpm run pbo_demo
fpm run --example basic_pbo
fpm run --example custom_metric
```

## Basic use

```fortran
use pbo, only : dp, pbo_result, compute_pbo, sharpe_ratio

type(pbo_result) :: result
real(dp) :: returns(1000, 50)

! Fill returns(time, strategy).
call compute_pbo(returns, 8, sharpe_ratio, result, threshold=0.0_dp)
if (.not. result%success) error stop result%message

print *, result%phi
print *, result%below_threshold
print *, result%degradation_slope
```

The performance argument is a procedure callback:

```fortran
subroutine metric(data, values)
  use pbo, only : dp
  real(dp), intent(in) :: data(:,:)
  real(dp), intent(out) :: values(:)
end subroutine metric
```

`data` contains one CSCV half-sample with observations in rows and candidate
strategies in columns. The callback returns one score per strategy. Built-in
callbacks are supplied for column means, column sums, zero-risk-free Sharpe
ratios, and zero-threshold Omega ratios. Parameterized variants are also
provided for direct use or from a user callback.

## Main result fields

- `combos`: subset combinations used by CSCV.
- `performance_is`, `performance_oos`: scores for every strategy and case.
- `selected_is`, `selected_oos`: first maximizing strategy indices.
- `oos_rank`: average OOS rank of the strategy selected in-sample.
- `omega_bar`: normalized OOS rank.
- `lambda`: rank logit, with positive infinity replaced by `inf_sub`.
- `phi`: probability of backtest overfitting.
- `selected_pairs(:,1:2)`: selected IS and corresponding OOS scores.
- `below_threshold`: probability that selected OOS performance is below the
  requested threshold.
- `degradation_intercept`, `degradation_slope`, `degradation_r2`: the directly
  interpretable OOS-on-IS degradation regression.

For compatibility, `slope`, `intercept`, and `adjusted_r2` reproduce the
original R object's regression calculation and naming, including its reversed
regression direction and swapped coefficient labels. See `PORTING_NOTES.md`.

## Plot-data replacements

The R package's lattice and grid graphics are not compiled. Their numerical
inputs remain available:

- `selection_frequencies` returns strategies sorted by decreasing selection
  count.
- `dominance_curve` returns the selected and all-strategy empirical CDFs, the
  upstream `SD2` pointwise difference, and an additional integrated difference.
- CSCV-case, pairs, rank, degradation, and lambda histogram data are already
  present in `pbo_result`.

## Licensing

The original package is MIT-licensed. The complete MIT license is in
`LICENSE`. Original source and supplied archive are retained under `original`
and `provenance`. See `NOTICE` and the checksum manifests.
