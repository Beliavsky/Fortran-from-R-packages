# Translation coverage

| Upstream functionality | Fortran status | Fortran entry point |
|---|---:|---|
| `fit_mfgarch` | Implemented | `fit_mfgarch` |
| Gaussian GARCH-MIDAS likelihood | Implemented | `log_likelihood`, `likelihood_contributions` |
| Simple GARCH with `K = 0` | Implemented | `fit_mfgarch` |
| GJR asymmetry | Implemented | `mfgarch_model%asymmetric` |
| Restricted beta weights | Implemented | `beta_weights`, `w1 = 1` |
| Unrestricted beta weights | Implemented | `unrestricted_weights = .true.` |
| One low-frequency covariate | Implemented | `covariate`, `period` |
| Two low-frequency covariates | Implemented | `has_second`, `covariate_two`, `period_two` |
| Irregular frequencies | Implemented | period-index arrays |
| Robust/HAC-style sandwich covariance | Implemented | `robust_covariance` |
| OPG covariance | Implemented | `opg_covariance` |
| BIC | Implemented | `mfgarch_fit_result%bic` |
| Variance ratio | Implemented | `variance_ratio` |
| Tau forecast | Implemented | `forecast_tau` |
| `predict.mfGARCH` numerical calculation | Implemented | `predict_variance` |
| `simulate_mfgarch` | Implemented | `simulate_mfgarch` |
| Student-t simulation | Implemented | optional `student_t_df` |
| `simulate_mfgarch_rv_dependent` | Implemented | `simulate_mfgarch_rv_dependent` |
| `simulate_mfgarch_diffusion` | Implemented | `simulate_mfgarch_diffusion` |
| C++ `calculate_g` | Implemented | `calculate_g` |
| C++ `sum_tau` / `sum_tau_fcts` | Implemented | same names |
| C++ `calculate_h_andersen` / `calculate_p` | Implemented | same names |
| C++ `simulate_r` | Implemented | same name |
| Plotting | Omitted | non-computational |
| R S3 classes and printing | Replaced | Fortran types and `print_fit_summary` |
| R data-frame/date validation | Replaced | array dimension/index validation |
| Bundled data objects | Not exposed as Fortran data | retained in `original/` |
