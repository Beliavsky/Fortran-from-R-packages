# RSDC: Regime-Switching Dynamic Correlation Models in R

**RSDC** is an R package for modeling time-varying asset correlations through a regime-switching framework. It supports flexible multivariate correlation dynamics with either fixed or time-varying transition probabilities (TVTP) driven by exogenous covariates.

This is particularly useful in financial applications where asset return correlations exhibit structural breaks or evolve under different market regimes.

## Features

- Three model families through one front end, `rsdc_estimate()`: constant
  correlations, fixed transition probabilities, and covariate-driven
  transitions (TVTP) with a logistic link for two regimes and a softmax for
  more
- Hamilton filter and smoother implemented in C++ (Rcpp/RcppArmadillo),
  validated against a pure-R reference to ~1e-8 in the test suite
- A global search that operates in a canonical partial-correlation space, in
  which every candidate is a valid positive-definite correlation matrix, so
  maximum likelihood stays feasible as the cross-section grows; plus
  data-driven warm starts (`rsdc_starts()`) and a multi-start replication
  diagnostic
- A full `"rsdc_fit"` S3 interface: `print`, `summary`, `coef`, `logLik`
  (so `AIC`/`BIC` work), `vcov`, `confint`, `predict`, `simulate`, `plot`,
  broom tidiers and `ggplot2::autoplot()`
- Standard errors from observed information, OPG, sandwich (QML) or a
  parametric bootstrap; uncertainty bands on the correlation path
- Forecasting and decoding: `rsdc_forecast()` (in-sample or strictly
  out-of-sample), `rsdc_forecast_ahead()` (multi-step), `rsdc_viterbi()`
- Portfolio construction from the implied covariances: `rsdc_minvar()`,
  `rsdc_maxdiv()`
- Shipped data: `greenbrown`, `mccc`, `ff5ind`

## Installation

```r
install.packages("RSDC")
```

Development version:

```r
devtools::install_github("ArdiaD/RSDC")
```

## Please cite the package in publications!

By using **RSDC**, you agree to the following:

1. Include a footnote or reference to the GitHub page:  
   [https://github.com/ArdiaD/RSDC](https://github.com/ArdiaD/RSDC)  
2. You assume all risk for using this software.

## References

Engle, R.F. (2002).
[Dynamic conditional correlation: A simple class of multivariate 
generalized autoregressive conditional 
heteroskedasticity models](https://doi.org/10.1198/073500102288618487).
*Journal of Business & Economic Statistics*, 20(3), 339–350.

Pelletier, D. (2006).
[Regime switching for dynamic correlations](https://doi.org/10.1016/j.jeconom.2005.01.013).
*Journal of Econometrics*, 131(1–2), 445–473.

Hamilton, J.D. (1989).
[A new approach to the economic analysis of nonstationary time 
series and the business cycle](https://doi.org/10.2307/1912559). 
*Econometrica*, 57(2), 357–384.  

## Companion article

A companion manuscript (in preparation for *The R Journal*) presents the
methodology, the estimation strategy, a Monte Carlo validation, and a strictly
out-of-sample study of five Fama–French industry portfolios with
regime-dependent minimum-variance and maximum-diversification allocations. It is
developed separately and is fully reproducible from data shipped with this
package. Use `citation("RSDC")` for the current reference.

For a short tour of the API, see the package vignette:
`browseVignettes("RSDC")`.
