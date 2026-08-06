# API map

| Upstream R/TMB entry point | Fortran entry point | Notes |
|---|---|---|
| `get_nll()` | `get_nll`, `laplace_nll`, `joint_nll` | Returns numerical objectives rather than a TMB AD object. |
| `estimate_parameters()` | `estimate_parameters` | Typed `sv_fit_result`; numerical Laplace and Nelder-Mead. |
| `sim_sv()` | `sim_sv` | All five implemented upstream models. |
| `simulate_parameters()` | `simulate_parameters` | Natural-scale parameter draws from transformed covariance. |
| `predict.stochvolTMB()` | `predict_sv` | Supports optional parameter uncertainty. |
| `summary.stochvolTMB_predict()` | `summarize_prediction` | Type-7 quantiles and optional means. |
| `summary.stochvolTMB()` | fields in `sv_fit_result` | Fixed, transformed, and latent estimates are directly accessible. |
| `logLik.stochvolTMB()` | `fit%log_likelihood` | Scalar field. |
| `residuals()` | `standardized_residuals`, `one_step_residuals` | Conditional PIT approximation at the latent mode. |
| `logit()` | `inv_logit_pm1` | Maps the real line to (-1,1). |
| TMB latent random effects | `latent_mode` | Sparse Newton mode and tridiagonal Hessian. |
| TMB `sdreport()` | `theta_cov`, `theta_se`, `param_se`, `h_se` | Finite-difference and Laplace uncertainty. |
| plotting functions | omitted | Non-computational. |
| `demo()` Shiny app | `app/demo_stochvoltmb.f90` | Terminal demonstration only. |
