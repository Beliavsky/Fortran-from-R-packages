# Array API

Use `use ftrading`; `dp` is re-exported from r_kinds. All computational routines
are pure. The scalar price transformations are elemental and accept conformable
arrays as well as scalars. Function outputs have the input series length except
trimmed rolling functions, trade lengths, and the matrices noted below.

## Public computational procedure signatures

Optional arguments and defaults are documented alongside each dummy in the
source, and summarized below. Arguments listed here preserve source order.

- `max_draw_down(x, drawdown, from, to)`
- `sharpe_ratio(x, risk_free, scale, rates) result(value)`
- `sterling_ratio(x) result(value)`
- `roll_var(x, n, trim, unbiased, na_rm) result(y)`
- `sma(x, n) result(y)`
- `ema_ta(x, lambda, startup) result(y)`
- `ewma(x, lambda, startup) result(y)`
- `trade_signals(positions, signals, indices)`
- `trade_lengths(times) result(lengths)`
- `hit_rate(returns, positions) result(rate)`
- `medprice_ta(high, low) result(y)`
- `typicalprice_ta(high, low, close) result(y)`
- `wclose_ta(high, low, close) result(y)`
- `vohl_ta(high, low) result(y)`
- `vor_ta(high, low) result(y)`
- `adoscillator_ta(open, high, low, close) result(y)`
- `mom_ta(x, lag) result(y)`
- `roc_ta(x, lag) result(y)`
- `bias_ta(x, lag) result(y)`
- `osc_ta(x, lag1, lag2) result(y)`
- `macd_ta(x, lag1, lag2) result(y)`
- `cds_ta(x, lag1, lag2, lag3) result(y)`
- `cdo_ta(x, lag1, lag2, lag3) result(y)`
- `fpk_ta(close, high, low, lag) result(y)`
- `wpr_ta(close, high, low, lag) result(y)`
- `fpd_ta(close, high, low, lag1, lag2) result(y)`
- `spd_ta(close, high, low, lag1, lag2, lag3) result(y)`
- `apd_ta(close, high, low, lag1, lag2, lag3, lag4) result(y)`
- `stochastic_ta(close, high, low, lag1, lag2, lag3, slow) result(y)`
- `rsi_ta(close, lag) result(y)`
- `accel_ta(x, n) result(y)`
- `adi_ta(high, low, close, volume) result(y)`
- `bollinger_ta(x, lag, n_sd) result(y)`
- `chaikino_ta(high, low, close, volume, lag1, lag2) result(y)`
- `chaikinv_ta(high, low, lag1, lag2) result(y)`
- `garmanklass_ta(open, high, low, close) result(y)`
- `nvi_ta(close, volume) result(y)`
- `pvi_ta(close, volume) result(y)`
- `obv_ta(close, volume) result(y)`
- `pvtrend_ta(close, volume) result(y)`
- `williamsad_ta(high, low, close) result(y)`
- `williamsr_ta(high, low, close, lag) result(y)`
- `daily_ta(ohlc, indicator, lags, select_column) result(y)`
- `roll_fun(x, n, operation, trim, na_rm)` accepts named `mean`, `sum`,
  `min`, `max`, `var`, `sd`, or `median` reductions.
- `roll_fun(x, n, fun, trim, na_rm)` accepts a procedure pointer conforming to
  the public pure `rolling_operation(window)` interface. Declare
  `procedure(rolling_operation), pointer :: reducer`, associate it with a pure
  scalar function, and pass the pointer. This permits explicit `intent(in)`
  on every procedure dummy. Extra callback arguments can be represented through
  a pure host-associated immutable configuration; R ellipsis is not reproduced.

## Shapes, defaults, and status conventions

- `max_draw_down` is a subroutine returning scalar magnitude plus allocatable
  `from(:)` / `to(:)` one-based index arrays. Every tied maximum is retained.
- `sharpe_ratio` defaults to scalar risk_free=0, scale=sqrt(250). Optional
  `rates(:)` must have size(x)-1 and replaces the scalar; no recycling occurs.
- `roll_fun` uses required positive n, trim=true, na_rm=false. Returned length
  is m-n+1 when trimmed, m otherwise, where m is the globally retained count.
  Invalid windows (n<1 or n>m) return size zero. Unknown named operations return
  NaNs. `roll_var` defaults n=9 and unbiased=true; `sma` defaults n=5.
- `ema_ta` defaults lambda=0.1, startup=0. `ewma` requires lambda. Decay >=1
  means period; decay=0 means full mean; startup=0 means floor(2/decay).
  Invalid or oversized startup returns same-shape NaNs.
- bias/ROC default lag=5; momentum=25; oscillator=(25,65); MACD=(12,26);
  CDS/CDO=(12,26,9); FPK=(5); FPD=(5,3); SPD=(5,3,9); APD=(5,3,9,9);
  WPR=50; RSI=14; acceleration n=12; Bollinger=(20,2); Chaikin oscillator
  periods=(10,3); Chaikin volatility=(10,10); Williams R lag=20.
- `stochastic_ta` returns an n-by-2 matrix, columns K and D. It defaults
  periods=(5,3,5), slow=false. Slow mode smooths both raw columns.
- `bollinger_ta` returns n-by-3: upper, original x, lower.
- `daily_ta` takes n-by-p OHLC data (p>=4) in Open/High/Low/Close order, defaults
  indicator='ema', lags=[9], select_column=4. Original named indicator choices
  are supported; missing required lag counts/invalid selectors yield NaNs.
- `trade_signals` returns allocatable nonzero full flips and matching one-based
  indices. `trade_lengths` differences numeric timestamps in the caller's unit.
  `hit_rate` drops NaN indicators and rounds to three decimals, ties-to-even.
  These four internal helper translations are auxiliary, outside the 40-function
  coverage denominator.
- Price/volume vectors must have equal length. Invalid shapes return same-shape
  NaNs. Empty arrays are safe; empty benchmark ratios are NaN. Comparisons are
  exact where upstream compares exact ties. NaNs propagate unless an explicit
  documented NA-removal argument applies. Zero divisions yield IEEE NaN/Inf.

The procedure declarations contain meaningful FORD `!!` documentation for each
argument's role, defaults, and shape. See README for upstream quirks and edge
cases where safely defined Fortran results replace R interface exceptions.
