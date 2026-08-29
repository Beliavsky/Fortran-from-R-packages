# API mapping

| R package computation | Fortran API | Notes |
|---|---|---|
| `dtweedie` | `dtweedie` | Scalar density; same interpolation/series selection |
| `ptweedie` | `ptweedie` | Scalar CDF; Fourier inversion outside closed forms |
| `qtweedie` | `qtweedie` | Scalar quantile by bracketed CDF inversion |
| `rtweedie` | `rtweedie`, `rtweedie_vec_params` | Compound Poisson-gamma for 1<p<2; inverse-CDF for p>2 |
| `dtweedie_series` | `dtweedie_series` | Dunn-Smyth series, both 1<p<2 and p>2 |
| `ptweedie_series` | `ptweedie_series` | Compound-Poisson series for 1<p<2 |
| `dtweedie_inversion` | `dtweedie_inversion` | Dunn-Smyth Fourier inversion, diagnostics available |
| `ptweedie_inversion` | `ptweedie_inversion` | Dunn-Smyth Fourier inversion, diagnostics available |
| `dtweedie_saddle` | `dtweedie_saddle` | Nelder-Pregibon offset behavior retained |
| `tweedie_dev` | `tweedie_dev` | Unit deviance |
| `tweedie_lambda` | `tweedie_lambda` | Compound-Poisson rate for p<2 |
| `tweedie_convert` | `tweedie_convert` | Returns `tweedie_convert_result` |
| `tweedie_integrand` | `tweedie_integrand_values` | Numerical real/imaginary/integrand values; plotting omitted |
| `logLiktweedie` | `tweedie_loglik` | Numeric likelihood for supplied y/mu/phi/power |
| `tweedie_AIC`, `AICtweedie` | `tweedie_aic` | Low-level AIC/general k-penalty calculation |
| internal `dtweedie_dlogfdphi` | `dtweedie_dlogfdphi` | Analytic/stabilized phi derivative |
| internal `dtweedie_dldphi` | `dtweedie_dldphi` | -2 log-likelihood derivative convention retained |
| profile phi optimization | `tweedie_phi_mle` | Positive scalar likelihood optimization |
| `tweedie_profile` core | `tweedie_glm_fit`, `tweedie_profile_grid` | Numeric design-matrix profile grid; no R formula/S3 layer |

The umbrella module is simply `use tweedie`.
