# API coverage audit

## Computational basis

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

## Per-function scope notes

- **ADX**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **ALMA**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **ATR**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **BBands**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **CCI**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **CLV**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **CMF**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **CMO**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **CTI**: Rolling Spearman correlation with average ranks, slope sign and window-local missingness. Zero/nonfinite slope gives NaNs.
- **DEMA**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **DPO**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **DVI**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **DonchianChannel**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **EMA**: SMA-seeded EMA with Wilder or explicit ratio. n remains explicit; inference of n from ratio and R warnings omitted.
- **EMV**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **EVWMA**: Elastic volume-weighted recurrence with first seed equal to current price. Scalar volume is explicitly expanded, repairing upstream wrapper row-count failure.
- **GMMA**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **HMA**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **KST**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **MACD**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **MFI**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **OBV**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **PBands**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **ROC**: Padded continuous or discrete lagged returns; zero-lag is a safe extension and nonpositive log-price inputs yield NaNs.
- **RSI**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **SAR**: Initial gap is sample SD of high-low range. Upstream equality-based acceleration cap, including possible overshoot, retained.
- **SMA**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **SMI**: Single common moving-average selector supported; upstream unbounded low-series runMax behavior retained. Per-stage R maType lists omitted.
- **SNR**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **TDI**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **TR**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **TRIX**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **VHF**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **VWAP**: Distinct exported alias of VWMA, exposed as its own public procedure.
- **VWMA**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **WMA**: Both fixed window weights and aligned series weights; oldest-to-newest order and upstream n-length priority retained.
- **WPR**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **ZLEMA**: Upstream fractional 1/ratio lag and SMA seed preserved; n remains explicit. Unsafe custom lags return NaNs.
- **ZigZag**: Percent/dollar, retracement, first/last extreme modes and linear interpolation; leading/trailing missing output retained.
- **adjRatios**: Split and dividend adjustment factors for caller-aligned arrays. Date merging and omission of missing-close rows are caller responsibilities.
- **aroon**: Upstream initial window, subsequent n+1 extrema window, latest-tie selection and leading-NA quirks retained.
- **chaikinAD**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **chaikinVolatility**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **growth**: Compounded ROC times signals. The default padded result propagates leading NaNs as upstream. Unpadded output is shortened without R vector recycling; signals must have matching length.
- **keltnerChannels**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **lags**: R embed ordering for time-by-variable matrices; incomplete leading rows omitted. Callers reshape a vector into one column.
- **momentum**: Padded lagged difference; zero lag yields zeros rather than R diff error. Oversized lag retains input-length NaNs.
- **naCheck**: Vector leading-NaN count and begin index plus validity flag; nonNA indices are the contiguous range n_missing+1:size(x). Matrix row checking is caller-managed; errors are flags.
- **rollSFM**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runCor**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runCov**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runMAD**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runMax**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runMean**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runMedian**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runMin**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runPercentRank**: Upstream asymmetric 1e-8 tie rule and cumulative initialization, including n=1 leading-missing quirk, retained.
- **runSD**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runSum**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **runVar**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted.
- **stoch**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **ultimateOscillator**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
- **volatility**: All six estimators, annualization and zero-mean/simple-return variants. Documented Yang-Zhang alpha/k overrides repaired; these throw object-not-found in upstream 0.24.4.
- **wilderSum**: Upstream seed x[n] plus (n-1)/n times the preceding sum, including n=1 leading-missing propagation, retained.
- **williamsAD**: Numerical array operation translated. R classes, indices, naming, warnings, and dynamic argument dispatch are omitted. Common maType selectors supported where applicable; R callback/list interfaces omitted.
