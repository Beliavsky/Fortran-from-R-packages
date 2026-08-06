# Porting notes

## Design choices

1. R reference classes were flattened into `trade_t`, a tagged modern Fortran
   derived type with type-bound procedures. This preserves computations while
   avoiding a large and brittle multiple-inheritance emulation.
2. R data-frame outputs became derived-type arrays. Lottery P&L output is held
   in `lottery_pnl_result_t`.
3. The `reticulate`/`pykalman` dependency was removed. `dynamic_beta` implements
   the same two-state random-walk regression model using a native Kalman filter,
   Rauch-Tung-Striebel smoother, and EM updates for the initial covariance,
   transition covariance, and observation variance.
4. The 139,838,160 Euro-lottery combinations are exposed through an iterator
   rather than allocated as one huge matrix.
5. Upstream plotting calls were deliberately omitted.

## Intentional correctness fixes

The upstream R source contains several apparent defects or unsafe edge cases.
The Fortran port applies these corrections and documents them rather than
silently reproducing them:

- `CrossSampleEntropy` now uses its `r` argument; the R implementation hardcodes
  `0.2` inside both matching loops.
- the UK Thunderball reader parses the actual textual draw-date column rather
  than referring to a non-existent `Date` column;
- the specific-number roulette simulation honors `targeted_number`; the R code
  always tests against zero;
- betting result arrays have exactly `simulations_number` entries and min/max
  statistics use the visited capital path rather than zero-filled unused cells;
- the bond split creates a `CDS`-equivalent `trade_t`; the R source calls an
  undefined `CreditSingle` constructor;
- CSA threshold processing avoids the final out-of-bounds index present in the
  R loop for some vector lengths;
- entropy routines validate input length and constant series rather than
  producing accidental indexing failures;
- lottery payout scales are calculated directly without allocating unused
  temporary payout vectors.

## Preserved upstream behavior that may be surprising

`martingale_strategy_repetitions` preserves the R algorithm's search for a run
of equal Boolean outcomes. Despite the R documentation saying "failed" trades,
the original algorithm counts runs of either successes or failures.

The roulette even-money strategies preserve the upstream rule that compares the
parity of consecutive spins. This is not the same as repeatedly betting a fixed
red/black selection.

## Numerical behavior

- sample variance and standard deviation use the `n - 1` denominator, matching
  R's `var` and `sd`;
- `quantile_type7` matches R's default type-7 quantile;
- the normal CDF uses the intrinsic complementary error function;
- covariance matrices are symmetrized and given small diagonal floors during
  Kalman EM iterations to avoid numerical singularity.

## Tests

`test/test_trading.f90` checks dependence metrics, type-7 quantiles, curve
interpolation, climate metrics, trade methods, digital-option CSV expansion,
curve/CSA/collateral/hash-table/track-record readers, lottery loading and payout
logic, combination iteration, betting simulations, and synthetic dynamic-beta
recovery.
