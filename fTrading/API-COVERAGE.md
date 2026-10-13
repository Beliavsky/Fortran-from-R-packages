# API coverage audit

Basis: distinct exported computational R functions, excluding plotting,
presentation, R startup hooks, and internal helpers. Package status substantial.
40 of 40 (100%) are mapped. This measures operations, not complete R compatibility.

Installed namespace has 45 exports: 40 user-facing computational functions,
four dot-prefixed internal helper leaks, and ohlcPlot. The 46 source-defined
functions also include the .onAttach hook. The manifest has exactly 40 mappings.

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

## Per-function compatibility notes

- `maxDrawDown` (upstream/R/BenchmarkAnalysis.R): All tied troughs and latest preceding peaks are returned as one-based indices. Empty or nonfinite levels yield NaN and empty indices; timeSeries timestamps omitted.
- `sharpeRatio` (upstream/R/BenchmarkAnalysis.R): Differences of levels, sample standard deviation, default sqrt(250) annualization; scalar or exactly aligned risk-free rates, without R recycling. Nonfinite levels rejected.
- `sterlingRatio` (upstream/R/BenchmarkAnalysis.R): Absolute-level drawdown denominator retained; zero denominators yield IEEE NaN or infinity; timeSeries metadata omitted.
- `rollFun` (upstream/R/rollFun.R): Pure scalar procedure-pointer callback or named scalar reduction, optional trim and global NaN removal. R ellipsis, vector-valued callbacks, class and index handling omitted.
- `rollVar` (upstream/R/rollFun.R): Sample variance or its (n-1)/n rescaling; global NaN removal and trim preserved. Invalid windows return an empty vector.
- `emaTA` (upstream/R/TechnicalAnalysis.R): First result is startup-window mean; periods map to decay 2/(period+1), default seed floor(2/decay). Zero decay is full-series mean. Invalid/oversized seed gives NaNs.
- `biasTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `medpriceTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `typicalpriceTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `wcloseTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `rocTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `oscTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `momTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `macdTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `cdsTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `cdoTA` (upstream/R/TechnicalAnalysis.R): Preserves upstream quirk: signal MACD uses default 12/26 periods even when the initial MACD has custom lag1/lag2.
- `vohlTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `vorTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `stochasticTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `fpkTA` (upstream/R/TechnicalAnalysis.R): Ratio in [0,1] for valid OHLC. Full window contains lag+1 observations, startup uses all available history.
- `fpdTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `spdTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `apdTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `wprTA` (upstream/R/TechnicalAnalysis.R): Preserves upstream definition exactly equal to fpkTA, rather than conventional Williams percent R; default lag is 50.
- `rsiTA` (upstream/R/TechnicalAnalysis.R): Unsmoothed sum-of-gains / sum-of-absolute-changes, width lag+1, scaled 0 to 1; first output copies second; flat window is NaN.
- `accelTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `adiTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `adoscillatorTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `bollingerTA` (upstream/R/TechnicalAnalysis.R): Columns upper, original x, lower; sample-standard-deviation bands and complete-window startup NaNs retained.
- `chaikinoTA` (upstream/R/TechnicalAnalysis.R): Array numerical operation; R classes, dates, labels, warnings and recycling omitted.
- `chaikinvTA` (upstream/R/TechnicalAnalysis.R): Numerical timeSeries path translated; vector-path upstream undefined rt bug repaired by using the computed range EMA.
- `garmanklassTA` (upstream/R/TechnicalAnalysis.R): Preserves unusual raw-price intraday differences with log overnight difference; not conventional scale-invariant Garman-Klass volatility.
- `nviTA` (upstream/R/TechnicalAnalysis.R): Cumulative selected ROC/100, where ROC uses current close denominator and falling volume selects observations. Initial value zero.
- `obvTA` (upstream/R/TechnicalAnalysis.R): Initial value zero and exact price ties contribute zero, unlike TTR initial-volume and near-tie conventions.
- `pviTA` (upstream/R/TechnicalAnalysis.R): Cumulative selected ROC/100, where ROC uses current close denominator and rising volume selects observations. Initial value zero.
- `pvtrendTA` (upstream/R/TechnicalAnalysis.R): Cumulative volume times percentage previous-close returns, initial zero; preserves upstream factor 100.
- `williamsadTA` (upstream/R/TechnicalAnalysis.R): Cumulative close-low on rising close, minus high-close on falling close; vector-path undefined x label bug omitted.
- `williamsrTA` (upstream/R/TechnicalAnalysis.R): Complete lag-width high/low windows and percentage distance below high; differs wprTA. Vector-path undefined x label bug omitted.
- `SMA` (upstream/R/TechnicalAnalysis.R): Trimmed complete-window arithmetic means; missing values affect their windows only; invalid windows return an empty vector.
- `EWMA` (upstream/R/TechnicalAnalysis.R): Same startup-seeded EMA operation as emaTA, with required decay/period argument. R class and index handling omitted.

## Exclusions

- `ohlcPlot`: plotting-only function.
- `.dailyTA`: namespace-pattern-exposed, undocumented internal helper.
- `.tradeSignals`: namespace-pattern-exposed, undocumented internal helper.
- `.tradeLengths`: namespace-pattern-exposed, undocumented internal helper.
- `.hitRate`: namespace-pattern-exposed, undocumented internal helper.
- `.onAttach`: package attach hook.
