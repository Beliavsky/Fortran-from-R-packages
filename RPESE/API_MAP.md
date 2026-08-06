# API map

| RPESE R function | Modern Fortran API | Notes |
|---|---|---|
| `EstimatorSE` | `estimate_se`, `estimate_se_matrix` | Generic vector and column-wise matrix dispatch. |
| `Mean.SE` | `mean_se` | All six SE methods. |
| `robMean.SE` | `robmean_se` | Robust cleaning is disabled, matching upstream intent. |
| `SD.SE` | `sd_se` | Sample standard deviation point estimate. |
| `SemiSD.SE` | `semisd_se` | Population-denominator downside dispersion. |
| `VaR.SE` | `var_se` | Wrapper accepts confidence and converts to lower-tail probability. |
| `ES.SE` | `es_se` | Wrapper accepts confidence and converts to lower-tail probability. |
| `SR.SE` | `sr_se` | Risk-free rate supported. |
| `SoR.SE` | `sor_se` | Mean or constant threshold. |
| `DSR.SE` | `dsr_se` | Source-compatible and corrected risk-free behavior. |
| `ESratio.SE` | `esratio_se` | Lower-tail probability and risk-free rate. |
| `VaRratio.SE` | `varratio_se` | Lower-tail probability and risk-free rate. |
| `RachevRatio.SE` | `rachevratio_se` | Separate lower and upper tail probabilities. |
| `LPM.SE` | `lpm_se` | Positive integer order; R interface documents 1 or 2. |
| `OmegaRatio.SE` | `omegaratio_se` | Threshold supported. |
| `SE.IF.iid` | `estimate_se(..., se_if_iid, ...)` | `sqrt(mean(IF^2)/n)`. |
| `SE.IF.cor` | `estimate_se(..., se_if_cor, ...)` | Periodogram elastic-net long-run variance. |
| prewhitened `SE.IF.cor` | `estimate_se(..., se_if_cor_pw, ...)` | AR(1) prewhitening and long-run correction. |
| `IFcorAdapt` | `estimate_se(..., se_if_cor_adapt, ...)` | Adaptive blend. |
| `SE.BOOT.iid` | `estimate_se(..., se_boot_iid, ...)` | Deterministic-seed iid bootstrap. |
| `SE.BOOT.cor` | `estimate_se(..., se_boot_cor, ...)` | Circular fixed-block bootstrap. |
| `fit.periodogram` | `fit_periodogram` | All, decimate, truncate, and optional two-sided output. |
| `SE.GLMEN` | `spectral_variance` | Exponential or Gamma RPEGLMEN fit. |
| `Add_Correlations` | fields in `se_result` | Return, raw IF, and prewhitened IF lag-one correlations. |
| `printSE` | omitted | Presentation-only R formatting. |
