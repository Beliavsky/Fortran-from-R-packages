# API map

| Original R/C routine | Modern Fortran counterpart | Status |
|---|---|---|
| `glag` | `lag_vector`, `lag_matrix`, generic `lag_values` | Implemented and tested |
| `gdiff` | `diff_vector`, `diff_matrix`, generic `diff_values` | Implemented and tested |
| `lgarchSim` / `LGARCHSIM` | `lgarch_simulate` | Implemented and tested, including arbitrary simulation orders |
| `lgarchRecursion1` / `ARMARECURSION1` | `lgarch_arma_recursion` | Implemented and tested, including zero imputation |
| `lgarchObjective` | `lgarch_objective` | LS, Gaussian QML, and CEX2 implemented and tested |
| `lgarch` | `fit_lgarch` | Implemented and tested for the original order restrictions |
| `coef.lgarch` | `lgarch_fit_result%arma_par`, `%lgarch_par` | Numerical result available; S3 method excluded |
| `fitted.lgarch` | `%fitted_sd`, `%log_sigma2` | Implemented and tested |
| `logLik.lgarch` | `%loglik_model`, `%loglik_arma`, `%objective_arma` | Implemented and tested |
| `residuals.lgarch` | `%residuals`, `%arma_residuals` | Implemented and tested |
| `rss` | `%rss` | Implemented and tested |
| `vcov.lgarch` | `%vcov_arma`, `%vcov_lgarch`, `%hessian_arma` | Implemented and tested with documented NaN blocks |
| `rmnorm` | `rmnorm` | Implemented and moment-tested |
| `mlgarchSim` | `mlgarch_simulate` | Implemented and tested |
| `mlgarchRecursion1` / `VARMARECURSION1` | `mlgarch_varma_recursion` | Implemented and tested |
| `mlgarchObjective` | `mlgarch_objective` | Implemented and checked against a direct likelihood calculation |
| `mlgarch` | `fit_mlgarch` | Implemented and tested for order zero or one |
| `coef.mlgarch` | `mlgarch_fit_result%varma_par`, `%mlgarch_par` | Numerical result available; S3 method excluded |
| `fitted.mlgarch` | `%fitted_sd`, `%log_sigma2` | Implemented and tested |
| `logLik.mlgarch` | `%objective_varma`, `%loglik_model` | Implemented and tested |
| `residuals.mlgarch` | `%residuals`, `%varma_residuals` | Implemented and tested |
| `vcov.mlgarch` | `%vcov_varma`, `%vcov_mlgarch`, `%hessian_varma` | Implemented and tested with documented NaN blocks |
| `print.*`, `summary.*` | None | R presentation/class infrastructure, excluded |
| `zoo` conversion and date-window behavior | Plain arrays and CSV date labels | R infrastructure, excluded |
| `c.code` switches | Not applicable | All recursions are compiled Fortran |
| Commented `rsep`, `rsst`, asymmetry, and prediction placeholders | None | Not implemented/exported by original package; not claimed |
