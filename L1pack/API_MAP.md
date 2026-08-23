# API map

| R / upstream entry | Fortran entry | Status |
|---|---|---|
| `dlaplace` | `dlaplace`, `dlaplace_vec` | translated |
| `plaplace` | `plaplace`, `plaplace_vec` | translated |
| `qlaplace` | `qlaplace`, `qlaplace_vec` | translated |
| `rlaplace` | `rlaplace`, `rlaplace_vec` | translated |
| `dmLaplace` | `dmlaplace`, `log_dmlaplace` | translated |
| `rmLaplace` | `rmlaplace` | translated |
| `l1fit` | `l1fit` | translated; Algorithm 478 |
| `lad.fit.BR` | `lad_fit_br` | translated |
| `lad.fit.EM` | `lad_fit_em` | translated |
| `lad.fit` | `lad_fit` | translated dispatcher |
| `lad` | `lad_fit` | numerical layer translated; R formula wrapper omitted |
| `vcov.lad` | `vcov_lad` | translated |
| `confint.lad` | `confint_lad` | translated |
| `predict.lad` | `predict_lad` | translated array API |
| `simulate.lad` | `simulate_lad` | translated |
| quantile residuals | `lad_quantile_residuals` | translated |
| `deviance.lad` | `lad_deviance` | translated |
| `logLik.lad` | `lad_result%loglik` | translated |
| `LaplaceFit` | `laplace_fit` | translated |
| internal equal-means fitter | `laplace_fit_equal` | translated |
| `spatial.median` | `spatial_median_fit` | translated |
| `l1ccc` | `l1ccc` | translated |
| Laplace CCC helper | `laplace_rho1` | translated |
| Gaussian CCC helper | `gaussian_rho1` | translated |
| U-statistic CCC helper | `ustat_rho1` | translated |
| bootstrap CCC variance | `l1ccc_bootstrap` | translated |
| `WH.Laplace` | `wh_laplace` | translated |
| `envelope.Laplace` | `envelope_laplace` | numerical envelope translated; plotting omitted |
| print/summary/plot S3 methods | — | R presentation layer omitted |
| formula/model-frame processing | — | R-specific layer omitted |
| bundled `ereturns` data | upstream under `orig/` only | data, not code |
