# Public Fortran API

Use `ttr` to import the complete library. Source dummy-argument FORD comments define shapes, units, defaults and valid ranges. All computation is pure. Optional arguments are specified in source; mandatory windows do not inherit R defaults.

## ttr_base

- `function nan_value()`
- `subroutine na_check(x, n, n_missing, begin, valid)`
- `integer function first_valid(x)`
- `function prices(x)`
- `subroutine high_low_close(x, h, l, c)`
- `function shift(x, n)`
- `function run_sum(x, n, cumulative)`
- `function run_mean(x, n, cumulative)`
- `function run_min(x, n, cumulative)`
- `function run_max(x, n, cumulative)`
- `function run_median(x, n, cumulative, non_unique)`
- `function run_cov(x, y, n, sample, cumulative)`
- `function run_var(x, n, sample, cumulative, y)`
- `function run_sd(x, n, sample, cumulative)`
- `function run_cor(x, y, n, sample, cumulative)`
- `function run_mad(x, n, center, median, constant, non_unique, cumulative)`
- `function run_percent_rank(x, n, cumulative, exact_multiplier)`
- `function wilder_sum(x, n)`
- `function sma(x, n)`
- `function ema(x, n, wilder, ratio)`
- `function dema(x, n, v, wilder, ratio)`
- `function wma(x, n, weights)`
- `function vwma(x, volume, n)`
- `function vwap(x, volume, n)`
- `function evwma(x, volume, n)`
- `function zlema(x, n, ratio)`
- `function hma(x, n)`
- `function alma(x, n, offset, sigma)`
- `function moving_average(x, n, method, default)`
- `function momentum(x, n)`
- `function roc(x, n, discrete)`
- `function lags(x, n)`
- `function growth(x, signals, discrete, n, na_pad)`

## ttr_indicators

- `function tr(hlc)`
- `function atr(hlc, n, ma_type)`
- `function adx(hlc, n, ma_type)`
- `function bbands(hlc, n, sd, ma_type)`
- `function cci(hlc, n, constant, ma_type)`
- `function clv(hlc)`
- `function cmf(hlc, volume, n)`
- `function cmo(x, n)`
- `function rsi(x, n, ma_type)`
- `function dpo(x, n, shift_by, percent, ma_type)`
- `function dvi(x, n, weights, smooth, magnitude, stretch, exact_multiplier)`
- `function donchian_channel(hl, n, include_lag)`
- `function emv(hl, volume, n, vol_divisor, ma_type)`
- `function gmma(x, periods, ma_type)`
- `function kst(x, n, n_roc, n_sig, weights, ma_type)`
- `function macd(x, n_fast, n_slow, n_sig, percent, ma_type)`
- `function mfi(hlc, volume, n)`
- `function obv(x, volume)`
- `function pbands(x, n, sd, fast_n, centered, long_average, ma_type)`
- `function snr(hlc, n, ma_type)`
- `function tdi(x, n, multiple)`
- `function trix(x, n, n_sig, percent, ma_type)`
- `function vhf(hlc, n)`
- `function wpr(hlc, n)`
- `function chaikin_ad(hlc, volume)`
- `function chaikin_volatility(hl, n, ma_type)`
- `function keltner_channels(hlc, n, atr_multiplier, ma_type)`
- `function stoch(hlc, n_fast_k, n_fast_d, n_slow_d, bounded, smooth, ma_type)`
- `function smi(hlc, n, n_fast, n_slow, n_sig, bounded, ma_type)`
- `function ultimate_oscillator(hlc, n, weights)`
- `function williams_ad(hlc)`

## ttr_advanced

- `function sar(hl, accel)`
- `function sar(high, low, accel)`
- `function zigzag(price, change, percent, retrace, last_extreme)`
- `function zigzag(hl, change, percent, retrace, last_extreme)`
- `function adj_ratios(splits, dividends, close)`
- `function aroon(price, n)`
- `function aroon(hl, n)`
- `function cti(price, n, slope)`
- `function volatility(price, n, calc, annual_periods, mean0, alpha, k, discrete)`
- `function volatility(ohlc, n, calc, annual_periods, mean0, alpha, k, discrete)`
- `function roll_sfm(asset, benchmark, n)`

## Matrix result columns (in order)

| Procedure | Columns |
| --- | --- |
| tr | tr, trueHigh, trueLow |
| atr | tr, atr, trueHigh, trueLow |
| adx | DIp, DIn, DX, ADX |
| bbands | dn, mavg, up, pctB |
| dvi | dvi.mag, dvi.str, dvi |
| donchian_channel | high, mid, low |
| emv | emv, maEMV |
| gmma | One column per supplied period, preserving order |
| kst | kst, signal |
| macd | macd, signal |
| pbands | dn, center, up |
| tdi | tdi, di |
| trix | TRIX, signal |
| keltner_channels | dn, mavg, up |
| stoch | fastK, fastD, slowD |
| smi | SMI, signal |
| aroon | aroonUp, aroonDn, oscillator |
| adj_ratios | Split, Div |
| roll_sfm | alpha, beta, r.squared |
| lags | All input columns at lag 0, then all at lag 1, through lag n |

All remaining mapped indicator functions return one vector. `sar`, `zigzag`, `aroon` and `volatility` are public generic interfaces: overloaded signatures above identify their accepted argument shapes. Helper routines such as `divide`, `prices`, `shift` and `moving_average` support composition but do not inflate R-function coverage. The listed public implementation utilities are additional Fortran composition helpers.
