# API map

## Direct numerical translations

| R routine | Fortran procedure |
|---|---|
| `dgev`, `.devd` | `gev_pdf`, `gev_logpdf` |
| `pgev`, `.pevd` | `gev_cdf` |
| `qgev`, `.qevd` | `gev_quantile` |
| `rgev`, `.revd`, `gevSim`, `gumbelSim` | `gev_random`, `gev_sample` |
| `gevMoments` | `gev_moments` |
| `dgpd`, `.depd` | `gpd_pdf`, `gpd_logpdf` |
| `pgpd`, `.pepd` | `gpd_cdf` |
| `qgpd`, `.qepd` | `gpd_quantile` |
| `rgpd`, `.repd`, `gpdSim` | `gpd_random`, `gpd_sample` |
| `gpdMoments` | `gpd_moments` |
| `.gumpwmFit` | `gumbel_pwm` |
| `.gevpwmFit` | `gev_pwm` |
| `.gummleFit`, `.gumLLH` | `fit_gumbel`, `gev_nll` |
| `.gevmleFit`, `.gevLLH`, `gevFit` | `fit_gev`, `gev_nll` |
| `.gpdpwmFit` | `gpd_pwm` |
| `.gpdmleFit`, `.gpdLLH`, `gpdFit` | `fit_gpd`, `gpd_nll` |
| `blockMaxima` (numeric blocks) | `block_maxima` |
| `findThreshold` | `find_threshold` |
| `pointProcess` | `point_process` |
| `deCluster` | `decluster` |
| `thetaSim` | `theta_simulate` |
| `blockTheta` | `block_theta` |
| `clusterTheta` | `cluster_theta` |
| `runTheta` | `run_theta` |
| `ferrosegersTheta` | `ferro_seg_theta` |
| `emdPlot` numerical output | `empirical_survival` |
| `qqparetoPlot` numerical output | `pareto_qq` |
| `mxfPlot`, `mePlot` numerical output | `mean_excess_curve` |
| `mrlPlot` numerical output | `mean_residual_life` |
| `recordsPlot` numerical output | `records_development` |
| `ssrecordsPlot` numerical output | `subsample_record_counts` |
| `msratioPlot` numerical output | `max_sum_ratios` |
| `sllnPlot` numerical output | `slln_path` |
| `lilPlot` numerical output | `lil_path` |
| `xacfPlot` numerical output | `exceedance_acf` |
| `hillPlot`, `shaparmHill` | `hill_estimator` |
| `shaparmPickands` | `pickands_estimator` |
| `shaparmDEHaan` | `dehaan_estimator` |
| `gpdTailPlot` numerical curve | `gpd_tail_curve` |
| `gpdQuantPlot`, `gpdShapePlot` numerical paths | `gpd_threshold_stability` |
| `gpdRiskMeasures`, `tailRisk` | `gpd_risk_measures` |
| `gpdQPlot`, `gpdSfallPlot`, `tailPlot` likelihood calculations | `gpd_profile_risk` |
| `gevrlevelPlot`, `.gevrlevelLLH` | `gev_return_level_profile` |
| `VaR`, `CVaR` sample modes | `sample_var`, `sample_cvar` |
| `.emaTA` | `ema_filter` |
| `.riskMetricsPlot` numerical filter | `riskmetrics_volatility` |
| `normMeanExcessFit` mean-excess formula | `normal_mean_excess` |

## Combined or plain-array equivalents

- `shaparmPlot`, `exindexPlot`, and `exindexesPlot` are plotting/orchestration
  wrappers. Their component estimators are callable separately.
- GEV fitting receives block maxima directly; callers may first invoke
  `block_maxima`.
- R fit objects are represented by `gev_fit_result`, `gpd_fit_result`,
  `risk_result`, `return_level_result`, and related plain derived types.

## Excluded dependency wrappers

- `ghMeanExcessFit`, `hypMeanExcessFit`, `nigMeanExcessFit`, and
  `ghtMeanExcessFit` call external `fBasics` fit/density/CDF routines.
- `.garch11MetricsPlot` calls external `fGarch::garchFit`.

## Excluded R infrastructure

Plotting and sliders, S3/S4 containers and methods, `timeSeries`/`timeDate`
metadata, package data attachment, formulas, and GUI behavior.
