# API mapping

| R / native entry point | Fortran counterpart | Notes |
|---|---|---|
| `nnet.default` | `nnet_fit` | Matrix API; formula/model-frame handling omitted. |
| `predict.nnet(type="raw")` | `nnet_predict` | Returns raw output matrix. |
| `predict.nnet(type="class")` | `nnet_predict_class` | Returns 1-based class indices. |
| `nnetHess` / `VR_nnHessian` | `nnet_hessian_exact` | Direct translation of analytic Hessian. |
| `VR_dfunc` | `nnet_objective_gradient` | Exact objective and analytic gradient. |
| `VR_nntest` | `nnet_predict_raw` | Low-level prediction with supplied weights. |
| `VR_set_net`, `norm.net`, `add.net` | `build_network` | Generates the standard topology used by exported `nnet`; arbitrary R list mutation is not reproduced. |
| `multinom` with factor response | `multinom_fit_labels` | Integer class labels `1..K`. |
| `multinom` with matrix response | `multinom_fit_counts` | Handles counts/proportions and censored category indicators. |
| `predict.multinom(type="probs")` | `multinom_predict_proba` | Probability matrix. |
| `predict.multinom(type="class")` | `multinom_predict_class` | 1-based class indices. |
| `multinomHess` | `multinom_information` | Equivalent baseline-category Fisher information. |
| `vcov.multinom` | `multinom_covariance` | LAPACK-SVD generalized inverse. |
| `logLik.multinom` | `multinom_loglik` | `edf`, deviance and AIC are fields of `multinom_model_t`. |
| `class.ind` | `class_ind` | One-hot matrix from integer labels. |
| `which.is.max` | `which_is_max` | Random tie breaking via Fortran RNG. |
| `VR_summ2` | `summarize_rows` | Lexicographic duplicate-X consolidation with Y summation. |

## R-specific code intentionally omitted

`nnet.formula`, model frames, contrasts, factor-level bookkeeping, S3
print/summary/coef wrappers, `add1`, `drop1`, `anova`, formula-based
`confint`, and other presentation/update machinery are not numerical kernels
and are not reproduced.  The quantities required for numerical post-processing
(weights, Hessians/information, covariance, deviance, AIC, log-likelihood,
rank, fitted values, and residuals) are retained.
