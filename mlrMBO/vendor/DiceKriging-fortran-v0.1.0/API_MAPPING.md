# API mapping

This table maps DiceKriging 1.6.1 computational functionality to the Fortran
API. R-only dispatch and presentation functions are listed separately.

| DiceKriging R API | Fortran counterpart | Notes |
|---|---|---|
| `km` | `km_fit`, `km_fit_cov`, `km_model` | R formula replaced by explicit trend matrix `f`. |
| `kmEstimate` | `km_estimate` | Native bounded multistart BFGS. |
| `computeAuxVariables` | `km_recompute` | Rebuilds Cholesky/whitened model state. |
| `kmNoNugget.init`, `km1Nugget.init`, `kmNuggets.init` | `km_fit_cov` fit-case initialization | Internalized rather than exposed as three R-style constructors. |
| `logLik.km` | `km_loglik` | Gaussian log likelihood for the fitted model. |
| `logLikFun` | `loglik_fun` | Profiled likelihood objective. |
| `logLikGrad` | `loglik_grad` | Analytic gradient for stationary kernels. |
| `leaveOneOut.km` | `leave_one_out` | SK/UK paths; exact slow fallback where required. |
| `leaveOneOutFun` | `leave_one_out_fun` | LOO squared-error criterion. |
| `leaveOneOutGrad` | `leave_one_out_grad` | Analytic LOO gradient. |
| `predict.km` / `predict` | `km_predict`, `km_prediction` | Mean, trend, SD, covariance, lower/upper 95% intervals. |
| `simulate` | `km_simulate` | Conditional or unconditional GP simulations. |
| `update` | `km_update` | `cov_reestimate=.true.` corresponds to covariance re-estimation. |
| `update(..., newX.alreadyExist=TRUE)` | `km_update_response` | Updates the most recently stored responses without rebuilding covariance. |
| `cv` | `cv_predict` | Fold labels are integer IDs. |
| `covTensorProduct` | `covariance_model` | `iso=.false.`, `scaling=.false.`. |
| `covIso` | `covariance_model` | `iso=.true.`. |
| `covScaling` | `covariance_model`, `scaling_axis` | Piecewise integrated scaling transform. |
| `covMatrix` | `covariance_matrix` | Supports nugget/noise. |
| `covMat1Mat2` | `covariance_cross` | Cross covariance, optionally matching nugget. |
| `covMatrixDerivative` | `covariance_derivative` | Range/shape/variance derivatives. |
| `covVector.dx` | `covariance_vector_dx`, `km_cov_vector_dx` | Spatial covariance derivative. |
| `covParametersBounds` | `covariance_bounds` | DiceKriging-style default bounds. |
| `covparam2vect` / `vect2covparam` | `get_cov_params` / `set_cov_params` | Flatten/unflatten covariance parameters. |
| `covStruct.create` | `covariance_model`, `km_fit` | Built-in covariance creation. |
| `scalingFun1d` | `scaling_fun1d`, `scaling_fun1d_dx` | 1-D scaling transformation and derivative. |
| `scalingFun` | `scaling_apply` | Multidimensional transformation. |
| `scalingGrad` | `scaling_grad1d` | Scaling-parameter gradient. |
| `trend.deltax` | `trend_gradient_*` | Helpers for common trend bases. |
| `trendMatrix.update` | explicit `newf` argument to `km_predict`/`km_update` | Formula parsing intentionally absent. |
| `SCAD` | `scad` | Elemental Fortran function. |
| `SCAD.derivative` | `scad_derivative` | Elemental derivative. |
| `branin` | `branin` | Native Fortran. |
| `camelback` | `camelback` | Native Fortran. |
| `goldsteinPrice` | `goldstein_price` | Native Fortran. |
| `hartman3` | `hartman3` | Native Fortran. |
| `hartman6` | `hartman6` | Native Fortran. |

## Trend helpers

- `trend_constant`: R formula `~ 1`
- `trend_linear`: intercept plus all first-order terms
- `trend_linear_interactions`: intercept, first-order terms, and all pairwise
  interactions (the `~ .^2` basis used by upstream tests)
- `trend_quadratic`: intercept, first-order terms, squares, and pairwise
  interactions

Matching `trend_gradient_*` routines return derivatives of these basis vectors.

## Intentionally omitted R/interface functionality

- `plot.km` and graphics imports: plotting code, explicitly outside the port.
- S4/S3 `show`, `summary`, `coef`, slot manipulation, class unions, and method
  dispatch: represented by public Fortran derived-type components and routines.
- `checkNames`, `checkNamesList`, `drop.response`, `kmData`: R
  data-frame/formula/name-management infrastructure.
- R `terms`, `model.matrix`, `reformulate`, and formula parsing: callers supply
  the numeric trend matrix directly.
- `covUser`: arbitrary R-function callback covariance. Built-in numerical
  kernels are translated; an R callback itself is not portable computational
  code.
- `rgenoud` and `foreach`/`doParallel`: external R dependencies, replaced by a
  serial native multistart bounded optimizer.
