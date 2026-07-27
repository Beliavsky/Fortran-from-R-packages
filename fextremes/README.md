# fExtremes modern Fortran

A modern Fortran translation of the computational parts of the R package
`fExtremes` 4032.84. The project uses plain arrays and derived result types;
plotting, sliders, S4 classes, formulas, and `timeSeries` metadata are not
reproduced.

## Implemented numerical functionality

### Probability distributions and simulation

- Generalized extreme-value (GEV) density, log density, CDF, quantile, moments,
  scalar RNG, and vector simulation
- Gumbel as the zero-shape GEV limit
- Generalized Pareto (GPD) density, log density, CDF, quantile, moments, scalar
  RNG, and vector simulation
- Stable support handling for positive, zero, and negative shape parameters

### Estimation and inference

- GEV and Gumbel probability-weighted-moment estimation
- GEV and Gumbel bounded maximum-likelihood estimation
- GPD probability-weighted-moment estimation
- GPD bounded maximum-likelihood estimation
- Observed numerical Hessians and covariance matrices
- Expected GPD covariance formula
- Standard errors and transformed exponential residuals
- GEV/GPD negative log-likelihood procedures

Bounded Nelder-Mead replaces R's `optim` orchestration. Exact optimizer paths
and endpoints are not claimed.

### Preprocessing and extremal dependence

- Numeric block maxima and source positions
- Threshold selection for a requested number of exceedances
- Point-process extraction with original indices
- Run declustering with cluster maxima and interval metadata
- Max-Frechet and paired-exponential extremal-index simulation
- Block, reciprocal-cluster-size, run, and Ferro-Segers extremal-index
  estimators

Calendrical monthly and quarterly block handling is excluded because it is
part of the R `timeDate`/`timeSeries` infrastructure. Callers can supply
numeric blocks or preprocess dates externally.

### Exploratory extreme-value calculations

- Empirical survival curves
- Exponential/GPD Pareto QQ coordinates
- Empirical mean-excess curves
- Mean residual life curves with normal confidence bands
- Record development and expected record counts
- Subsample record counts
- Maximum-to-sum ratio paths
- Strong-law and iterated-logarithm paths
- Exceedance heights, inter-exceedance distances, and their autocorrelations
- Pickands, Hill, and Dekkers-Einmahl-de Haan tail-index curves and summaries
- Normal-distribution mean excess

### Risk and return levels

- Sample lower- or upper-tail VaR and CVaR
- GPD tail survival curves
- GPD VaR and expected shortfall with delta-method standard errors
- GPD profile-likelihood intervals for VaR and expected shortfall
- Threshold-stability paths for shape, scale, modified scale, uncertainty, and
  a selected high quantile
- GEV return levels with delta-method intervals
- GEV profile-likelihood return-level intervals
- EMA filtering and RiskMetrics-style volatility

## Applications

```text
fit_csv FILE gev
fit_csv FILE gumbel
fit_csv FILE gpd [threshold]
```

The CSV reader accepts one numeric column, or a date/text column followed by a
numeric column, with comma or semicolon separators.

```text
./build/debug/demo_fextremes
./build/debug/fit_csv data/danishClaims.csv gpd 10
./build/debug/tail_analysis
```

## Build and test

```sh
make debug
make release
make check
```

The tested GNU Fortran builds use Fortran 2018, warnings as errors, and full
runtime checks in debug mode. LAPACK/BLAS are not required.

An `fpm.toml` file is included, but `fpm` was unavailable in the validation
environment and is not claimed as tested.

## Explicit omissions

- Plotting, interactive sliders, S3/S4 classes, print/summary methods, formulas,
  and R date/time-series metadata
- The GH, hyperbolic, NIG, and GH-Student mean-excess fitting wrappers. Those
  functions delegate distribution fitting and integration to `fBasics`; the
  original `fExtremes` package does not contain their numerical engines.
  A normal mean-excess implementation is included.
- The internal GARCH(1,1) volatility plotting wrapper, which delegates fitting
  to `fGarch`. EMA and RiskMetrics calculations are included.
- Exact R spline interpolation, `optim`, `uniroot`, random streams, and
  iteration-by-iteration equivalence

These omissions are dependency or R-infrastructure boundaries, not silently
claimed features.
