# R-to-Fortran API map

## Exported computational functions

| R function | Modern Fortran procedure | Notes |
|---|---|---|
| `B.specFGN` | `b_spec_fgn` | `k_approx <= 0` selects direct summation. |
| `CetaARIMA` | `ceta_arima` | Returns an allocatable covariance matrix. |
| `CetaFGN` | `ceta_fgn` | Intended for the one-parameter fGn model. |
| `ckARMA0` | `ck_arma0` | fARIMA(0,d,0) autocovariances. |
| `ckFGN0` | `ck_fgn0` | fGn autocovariances. |
| `FEXPest` | `fexp_estimate` | Returns `type(fexp_result)`. |
| `per` | `periodogram` | Includes the zero-frequency ordinate. |
| `.ffreq` | `fourier_frequencies` | Fourier frequencies excluding zero. |
| `Qeta` | `qeta` | Returns `type(qeta_result)` or only B. |
| `simARMA0` | `sim_arma0` | Gaussian simulation from fARIMA(0,d,0) covariance. |
| `simFGN0` | `sim_fgn0` | Gaussian simulation from fGn covariance. |
| `simFGN.fft` | `sim_fgn_fft` | Paxson spectral simulation; requires even n. |
| `simGauss` | `sim_gauss` | Circulant-embedding Gaussian simulation. |
| `specARIMA` | `spec_arima` | Returns `type(spectrum_result)`. |
| `specFGN` | `spec_fgn` | Returns `type(spectrum_result)`. |
| `WhittleEst` | `whittle_estimate` | Returns `type(whittle_result)`. |

## Result-object methods

R's S3 accessors are represented directly by derived-type fields:

| R method | Fortran equivalent |
|---|---|
| `coef.FEXP` | `fexp_result%coefficients(:,1)` |
| `coef.WhittleEst` | `whittle_result%eta` |
| `nobs.FEXP` | `fexp_result%n` |
| `nobs.WhittleEst` | `whittle_result%n` |
| `vcov.FEXP` | `fexp_result%covariance` |
| `vcov.WhittleEst` | `whittle_result%covariance` |
| `print.FEXP` | ordinary formatted output by the caller |
| `print.WhittleEst` | ordinary formatted output by the caller |

For `fexp_result%coefficients`, columns are:

1. estimate
2. standard error
3. t value
4. two-sided p value

## Deliberately omitted plotting code

The following graphics-only routines and methods were not translated:

- `llplot`
- `lxplot`
- `lines.FEXP`
- `lines.WhittleEst`
- `plot.FEXP`
- `plot.WhittleEst`

The numerical frequency, periodogram, and fitted-spectrum arrays remain available for plotting with any external tool.
