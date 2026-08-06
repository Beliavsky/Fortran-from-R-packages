# API map

| Upstream R/C++ routine | Modern Fortran routine | Notes |
|---|---|---|
| `glmnet_exp` / `fitGlmCv` | `glmnet_exp` / `fitGlmCv` | Exponential elastic-net CV |
| `fit.glmGammaNet` | `fit_glm_gamma_net` / `fit_glmGammaNet` | Selects Gamma or exponential model |
| `cv.glmGammaNet` | `fit_glm_gamma_net` / `cv_glmGammaNet` | Gamma shape MLE followed by CV |
| `glmGammaNet` | `glm_gamma_net_fixed` / `glmGammaNet` | Fixed shape and lambda |
| `fit_glm_gamma` | `fit_glm_gamma_mle` | Joint coefficients and log-shape BFGS |
| `FISTA` / `ProxGradDescent` | `fit_fixed_model` | Backtracking proximal solver |
| `prox_L1` | `prox_l1` | Optional unpenalized first coefficient |
| `prox_EN` | `prox_en` | Corrected and source-compatible forms |
| `regularizer_EN` | `regularizer_en` | L1/L2 elastic-net penalty |
| `mse` | `mse_value` | Half squared-error loss |
| `mse_grad` | `mse_gradient` | Squared-error gradient |
| `nll_gamma_glm` | `gamma_negative_log_likelihood` | Fixed-shape mean NLL |
| `nll_gamma_glm_grad` | `grad_gamma_negative_log_likelihood` | Analytic gradient |
| `nll_gamma_optim` | internal `gamma_joint_objective_gradient` | Uses log-shape parameterization |
| `ExpNegativeLogLikelihood` | `exp_negative_log_likelihood` | Source uses an unnormalized sum |
| `GradExpNegativeLogLikelihood` | `grad_exp_negative_log_likelihood` | Analytic gradient |
| `Predict` | `predict_mean` | Log-link inverse, `exp(A beta)` |
| `ComputeLambdaMax` | `compute_lambda_max` | Handles unpenalized intercept and ridge |
| `GenerateLambdaGrid` | `generate_lambda_grid` | Descending logarithmic grid |
| `FitGlmFixed` | `fit_regularization_path` or fixed-lambda routines | Upstream C++ routine is a zero-vector stub |
| `GenerateCvData` | internal balanced fold construction | Corrects upstream integer-index defects |
| `scalarMultiplication`, `addReals`, `MyClass` | omitted | Rcpp/Eigen interface demonstrations only |
