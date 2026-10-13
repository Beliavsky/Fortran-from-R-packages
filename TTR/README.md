# TTR: Technical Trading Rules in modern Fortran

A substantial computational translation of **TTR 0.24.4** by Joshua M. Ulrich and contributors. All maintained Fortran is free-form, pure, and uses the single shared `dp` kind. This is a numerical library, not investment advice or a trading recommendation.

## Build and run

Place `TTR/` alongside the repository's existing `rfortran-core/`:

```text
Fortran-from-R-packages/
  TTR/
  rfortran-core/
```

```sh
cd TTR
fpm build
fpm test
fpm run --example
```

Only FPM and a Fortran 2008-capable compiler such as gfortran are needed, including on Windows. No R installation, system BLAS/LAPACK, network access, or external library linker flags are needed to build or run. The sibling dependency is not bundled. Do not enable `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent flags: missing values are part of the API.

Import `use ttr` for the public API, or select `ttr_base`, `ttr_indicators`, and `ttr_advanced`. Public arguments are documented inline using FORD comments; [API.md](API.md) lists signatures and result column ordering.

```fortran
use ttr, only: dp, ema, rsi
real(dp) :: close(100), average(100), strength(100)
! Populate close chronologically before these calls.
average = ema(close, 10)
strength = rsi(close, 14)
```

## Array and numerical conventions

- Time is the first dimension, oldest observation first. Prices, volumes, returns and windows must already be aligned. Single-series operations take vectors. Indicator matrices use HLC = high/low/close, HL = high/low, OHLC = open/high/low/close. Functions accepting HLC or close take a one-column matrix for the close form.
- Results generally preserve input rows and use IEEE quiet NaN for warm-up/missing output. `lags` drops incomplete rows; `growth(...,na_pad=.false.)` returns a shortened allocatable vector. `adj_ratios` returns an allocatable two-column matrix.
- Window lengths are explicit for most routines; do not assume R defaults. Optional parameters retain their documented defaults. Lagged returns are always padded except the `growth` option. No vector recycling, R classes, timestamps, names, S3 reconstruction, lazy evaluation, warning delivery, or implicit date joins are implemented.
- Most rolling/recursive operations allow a leading missing prefix and reject interior missing data by returning all-NaN outputs. Invalid windows, incompatible shapes, unsupported selectors, or insufficient observations return NaNs rather than throwing R exceptions. `na_check` returns a validity flag and counts; its vector API can be applied to a caller-computed row-missingness mask. CTI is explicitly window-local and can resume after interior NaNs leave its window. Algebraic elementwise functions propagate missingness locally. R NA and computational NaN share one IEEE representation, so exceptional-data distinctions are not exactly reproducible.
- Moving-average `ma_type` selectors are SMA, EMA, Wilder, WMA, DEMA, HMA, ALMA and ZLEMA (lowercase accepted). A single selector applies to all stages; arbitrary R callbacks, per-stage lists and additional MA-specific ellipsis arguments are omitted. Apply the exposed primitives directly when a custom pipeline is needed.
- Rolling sums preserve upstream C seed/add/subtract order, important for rank-based indicators. Covariance, variance, median and sorting reuse compatible `rfortran-core` routines; small floating-point differences remain possible. EMA seeds use a mean of the first window. Wilder sum has the distinct upstream seed, not a plain initial sum.

## Intentional differences and retained quirks

- `EMA`, `DEMA` and `ZLEMA` require `n`, even when a decay ratio is supplied. No ratio-to-window inference or R warning is emitted.
- Scalar EVWMA volume is expanded explicitly; the R wrapper advertises this but fails its matrix row check. Unpadded `growth` avoids accidental R default-signal recycling and requires signals to match its output length. Zero-lag momentum and ROC are safe extensions; oversized lags return input-length NaNs.
- Yang–Zhang volatility supports the documented `alpha` and `k` overrides, repairing upstream 0.24.4 object-not-found errors. Those two expected results are independently evaluated using the documented formula and genuine TTR component operations.
- The upstream runPercentRank cumulative initialization/tie comparison, wilderSum `n=1` leading-NA propagation, OBV logarithmic-change tolerance and cumulative missing propagation, SMI unbounded low-series `runMax`, SAR equality-based acceleration cap, aroon initialization/subsequent `n+1` range, and ZigZag trailing missingness are deliberately retained.
- Adjustment ratios assume event rows are already aligned with prices; callers must merge dates and remove missing-close rows themselves. Matrix `naCheck` dispatch and generated R column names are outside the array API.
- Nonpositive prices are outside the log-return/volatility mathematical domain. Degenerate divisions use IEEE NaN/infinity, with CLV, WPR, MFI and stochastic special cases matching their documented numeric conventions. Floating-point exception-trapping builds are not supported for intentionally missing/degenerate data.

## Verification

`fpm test` runs 89 main genuine-R cases (12,941 numerical entries), 42 advanced cases (1,488 entries), a 512-row DVI rounding regression (1,536 expected entries), a direct naCheck fixture, and independent edge/algebraic properties. Coverage includes every mapped export, all six volatility estimators, stateful reversals, repeated extrema, tie rules, custom weights, cumulative windows, missing prefixes, initialization and matrix column ordering. The tests do not call R at runtime.

Reproduce stored fixtures with genuine TTR 0.24.4 installed:

```sh
Rscript test/reference.R
Rscript test/advanced_reference.R
Rscript test/dvi_reference.R
fpm test
```

The generators and expected results are committed; they are independent reference calculations against the unmodified R package, not outputs captured from this Fortran implementation. See [NOTICE.md](NOTICE.md) for source provenance and licensing.

## Translation coverage

**Status: substantial. 65 of 65 (100%)** exported computational R functions have meaningful mappings. This fraction measures mapped computational functions, not complete R compatibility. All mappings are substantial because interface and error behavior differ. The denominator contains the 67 explicit exports minus data-download functions `getYahooData` and `stockSymbols`. `naCheck` is included conservatively because it computes missing-count and index information. The distinct `VWAP` and `VWMA` exported names are both counted, including their alias relationship. There are no untranslated functions in this computational basis.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| ADX | ttr_indicators | adx | substantial |
| ALMA | ttr_base | alma | substantial |
| ATR | ttr_indicators | atr | substantial |
| BBands | ttr_indicators | bbands | substantial |
| CCI | ttr_indicators | cci | substantial |
| CLV | ttr_indicators | clv | substantial |
| CMF | ttr_indicators | cmf | substantial |
| CMO | ttr_indicators | cmo | substantial |
| CTI | ttr_advanced | cti | substantial |
| DEMA | ttr_base | dema | substantial |
| DPO | ttr_indicators | dpo | substantial |
| DVI | ttr_indicators | dvi | substantial |
| DonchianChannel | ttr_indicators | donchian_channel | substantial |
| EMA | ttr_base | ema | substantial |
| EMV | ttr_indicators | emv | substantial |
| EVWMA | ttr_base | evwma | substantial |
| GMMA | ttr_indicators | gmma | substantial |
| HMA | ttr_base | hma | substantial |
| KST | ttr_indicators | kst | substantial |
| MACD | ttr_indicators | macd | substantial |
| MFI | ttr_indicators | mfi | substantial |
| OBV | ttr_indicators | obv | substantial |
| PBands | ttr_indicators | pbands | substantial |
| ROC | ttr_base | roc | substantial |
| RSI | ttr_indicators | rsi | substantial |
| SAR | ttr_advanced | sar | substantial |
| SMA | ttr_base | sma | substantial |
| SMI | ttr_indicators | smi | substantial |
| SNR | ttr_indicators | snr | substantial |
| TDI | ttr_indicators | tdi | substantial |
| TR | ttr_indicators | tr | substantial |
| TRIX | ttr_indicators | trix | substantial |
| VHF | ttr_indicators | vhf | substantial |
| VWAP | ttr_base | vwap | substantial |
| VWMA | ttr_base | vwma | substantial |
| WMA | ttr_base | wma | substantial |
| WPR | ttr_indicators | wpr | substantial |
| ZLEMA | ttr_base | zlema | substantial |
| ZigZag | ttr_advanced | zigzag | substantial |
| adjRatios | ttr_advanced | adj_ratios | substantial |
| aroon | ttr_advanced | aroon | substantial |
| chaikinAD | ttr_indicators | chaikin_ad | substantial |
| chaikinVolatility | ttr_indicators | chaikin_volatility | substantial |
| growth | ttr_base | growth | substantial |
| keltnerChannels | ttr_indicators | keltner_channels | substantial |
| lags | ttr_base | lags | substantial |
| momentum | ttr_base | momentum | substantial |
| naCheck | ttr_base | na_check | substantial |
| rollSFM | ttr_advanced | roll_sfm | substantial |
| runCor | ttr_base | run_cor | substantial |
| runCov | ttr_base | run_cov | substantial |
| runMAD | ttr_base | run_mad | substantial |
| runMax | ttr_base | run_max | substantial |
| runMean | ttr_base | run_mean | substantial |
| runMedian | ttr_base | run_median | substantial |
| runMin | ttr_base | run_min | substantial |
| runPercentRank | ttr_base | run_percent_rank | substantial |
| runSD | ttr_base | run_sd | substantial |
| runSum | ttr_base | run_sum | substantial |
| runVar | ttr_base | run_var | substantial |
| stoch | ttr_indicators | stoch | substantial |
| ultimateOscillator | ttr_indicators | ultimate_oscillator | substantial |
| volatility | ttr_advanced | volatility | substantial |
| wilderSum | ttr_base | wilder_sum | substantial |
| williamsAD | ttr_indicators | williams_ad | substantial |
