# fTrading: modern Fortran / FPM

Numerical translation of CRAN fTrading 3042.79 by Diethelm Wuertz, Tobias Setz,
and Yohan Chalabi. The original package has GPL (>= 2) metadata; its R source
files contain GNU Library GPL version 2-or-later notices. Both are retained.
See [NOTICE.md](NOTICE.md), [LICENSE](LICENSE), and [LICENSES](LICENSES).

## Build, test, and example

Place this directory at the top level of Fortran-from-R-packages beside `TTR`
and `rfortran-core`, then run:

```
fpm build
fpm test
fpm run --example trading_example
```

All maintained Fortran is free form and uses `dp` from `r_kinds`. No separately
installed BLAS/LAPACK or system library link flags are required. The package
uses only portable Fortran 2008 facilities, so FPM with gfortran on Windows
needs the same sibling packages, not a Unix library installation. R is needed
only to regenerate fixtures, not to build, test, or run the package.

Import `use ftrading` for the public API. Every numerical routine is pure;
pointwise price transformations are pure elemental. See [API.md](API.md) for
signatures, shapes, defaults, and missing-value rules.

## Compatibility and deliberate details

- Inputs are chronological `real(dp)` arrays. Dates, `timeSeries` objects,
  plotting, S3/R object construction, names, automatic vector recycling,
  dynamic argument dispatch, and R warnings/errors are omitted.
- `max_draw_down` returns an absolute level decline, all tied troughs, and
  the latest preceding peak for each trough. Indices start at one. Sharpe and
  Sterling consume levels (often log prices), not a ready-made return vector.
- EMA is intentionally different from TTR's ordinary initialization. The first
  result equals the mean of the startup observations; recursion starts at the
  second observation. A period >=1 becomes decay 2/(period+1), and default
  startup is floor(2/decay). Zero decay produces the full-series mean.
- `fpk_ta` uses lag+1 observations and a partial startup window. `wpr_ta` is
  identical to `fpk_ta`, as the original code specifies. `rsi_ta` uses a
  lag+1 unsmoothed window, produces a 0-to-1 ratio, and copies result 2 to 1.
- `cdo_ta` preserves the original call: its signal line always uses default
  MACD periods 12 and 26 even with custom lag1/lag2.
- Bollinger columns are upper, original input, lower. The middle column is
  not the moving mean. Bands use sample standard deviation.
- `garmanklass_ta` preserves raw-price intraday differences plus a logarithmic
  overnight term. Its output is scale-dependent and should not be confused
  with the conventional log-price volatility estimator.
- NVI/PVI are cumulative selected current-close-denominator ROC divided by
  100. Price-volume trend uses previous-close percentage changes. OBV starts
  at zero and treats only exact equal prices as ties.
- `chaikinvTA` has an undefined `rt` on its R vector path; `williamsadTA` and
  `williamsrTA` reference undefined `x` when adding vector names. Their working
  timeSeries numerical paths are translated and exposed safely to arrays.
- NaNs represent R NA/NaN. Rolling windows recover after an interior NaN;
  `na_rm` removes NaNs globally before constructing windows. EMA propagates
  NaNs through its seed/recursion. Zero ratios return IEEE NaN or signed
  infinity. Do not compile with fast-math, Ofast, finite-math-only, or similar
  flags that discard NaN/infinity semantics.
- Invalid vector lengths/periods generally return same-shape NaNs; invalid
  rolling windows return an empty vector. Empty/nonfinite benchmark inputs
  return NaN (and empty drawdown indices). R errors, malformed negative-lag
  slicing, recycling, and accidental source exceptions are not reproduced.
  Arithmetic overflow may produce IEEE special values; rescaled mean
  and variance fallbacks handle finite constant/large-value windows when direct
  shared-kernel accumulation overflows.

## Dependency reuse audit

The existing top-level TTR and rfortran-core APIs were inspected before
implementation. `TTR = { path = "../TTR" }` and
`rfortran-core = { path = "../rfortran-core" }` are sibling dependencies.
No source from either dependency is copied into this package.

- TTR momentum is reused, adapting its leading NaNs to fTrading zero padding.
- TTR EMA is reused with an adapted first observation and one-value seed;
  a narrow nonfinite fallback preserves fTrading recursive propagation,
  which TTR's leading-NA-only precondition does not permit.
- TTR rolling mean/min/max kernels operate on complete individual windows,
  preserving fTrading window-local NA behavior. Core variance, standard
  deviation, mean, and median implement compatible scalar statistics.
- TTR ROC, RSI, stochastic, OBV, Bollinger, and volume indicators have different
  denominators, seed points, window widths, scaling, or tie behavior. Their
  mathematical definitions are not interchangeable; only compatible kernels
  are reused, with package-specific formulas preserved above.
- TTR's explicit IEEE-safe divide and NaN helpers are reused. No linear-algebra
  dependency is needed. R-only timeDate/timeSeries/fBasics dependencies are not
  reimplemented, vendored, or required by the Fortran API.

## Verification

`test_invariants` independently checks drawdown tie indices, sample-return
ratios, rolling arithmetic and variance, pure callback support, NaN recovery,
invalid lengths/windows, EMA recurrence, large magnitudes, RSI/stochastic
bounds, factor/denominator conventions, volume partitions, and auxiliary
trade operations. R reference fixtures are generated from the genuine
installed fTrading package and cover the mapped computational API, including
known interface-bug functions through their working timeSeries paths.
See `test/reference` for fixture generation and exact provenance. The example
is deterministic. Default, checked debug, release, example, and clean-extraction
validation are run separately by the batch verification workflow.

## Translation coverage

Package status: **substantial**. **40 of 40 (100%)** exported computational
R functions are mapped. This fraction measures mapped computational functions,
not complete R compatibility. All mapping statuses are substantial because
the R object/interface layer is intentionally omitted.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| maxDrawDown | ftrading_base | max_draw_down | substantial |
| sharpeRatio | ftrading_base | sharpe_ratio | substantial |
| sterlingRatio | ftrading_base | sterling_ratio | substantial |
| rollFun | ftrading_base | roll_fun | substantial |
| rollVar | ftrading_base | roll_var | substantial |
| emaTA | ftrading_base | ema_ta | substantial |
| biasTA | ftrading_indicators | bias_ta | substantial |
| medpriceTA | ftrading_indicators | medprice_ta | substantial |
| typicalpriceTA | ftrading_indicators | typicalprice_ta | substantial |
| wcloseTA | ftrading_indicators | wclose_ta | substantial |
| rocTA | ftrading_indicators | roc_ta | substantial |
| oscTA | ftrading_indicators | osc_ta | substantial |
| momTA | ftrading_indicators | mom_ta | substantial |
| macdTA | ftrading_indicators | macd_ta | substantial |
| cdsTA | ftrading_indicators | cds_ta | substantial |
| cdoTA | ftrading_indicators | cdo_ta | substantial |
| vohlTA | ftrading_indicators | vohl_ta | substantial |
| vorTA | ftrading_indicators | vor_ta | substantial |
| stochasticTA | ftrading_indicators | stochastic_ta | substantial |
| fpkTA | ftrading_indicators | fpk_ta | substantial |
| fpdTA | ftrading_indicators | fpd_ta | substantial |
| spdTA | ftrading_indicators | spd_ta | substantial |
| apdTA | ftrading_indicators | apd_ta | substantial |
| wprTA | ftrading_indicators | wpr_ta | substantial |
| rsiTA | ftrading_indicators | rsi_ta | substantial |
| accelTA | ftrading_indicators | accel_ta | substantial |
| adiTA | ftrading_indicators | adi_ta | substantial |
| adoscillatorTA | ftrading_indicators | adoscillator_ta | substantial |
| bollingerTA | ftrading_indicators | bollinger_ta | substantial |
| chaikinoTA | ftrading_indicators | chaikino_ta | substantial |
| chaikinvTA | ftrading_indicators | chaikinv_ta | substantial |
| garmanklassTA | ftrading_indicators | garmanklass_ta | substantial |
| nviTA | ftrading_indicators | nvi_ta | substantial |
| obvTA | ftrading_indicators | obv_ta | substantial |
| pviTA | ftrading_indicators | pvi_ta | substantial |
| pvtrendTA | ftrading_indicators | pvtrend_ta | substantial |
| williamsadTA | ftrading_indicators | williamsad_ta | substantial |
| williamsrTA | ftrading_indicators | williamsr_ta | substantial |
| SMA | ftrading_base | sma | substantial |
| EWMA | ftrading_base | ewma | substantial |

The broad namespace pattern also exposes `.dailyTA`, `.tradeSignals`,
`.tradeLengths`, and `.hitRate`. These undocumented internal helpers are
excluded from the denominator under the internal-helper rule, but their
numerical operations are provided as `daily_ta`, `trade_signals`,
`trade_lengths`, and `hit_rate`. `ohlcPlot` (plotting) and `.onAttach`
(startup presentation hook, not exported by the installed namespace) are
excluded. Imported functions and commented-out `rollMean`, `rollMin`, and
`rollMax` definitions are not counted. No user-facing computational export
is left untranslated. See [API-COVERAGE.md](API-COVERAGE.md) for the complete
inventory and provenance paths.
