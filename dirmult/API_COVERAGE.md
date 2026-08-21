# API coverage

| R routine | Fortran equivalent | Status |
|---|---|---|
| `dirmult` | `fit_dirmult` | translated |
| `dirmult.summary` | `summarize_dirmult` | translated |
| `weirMoM` | `weir_mom` | translated |
| `estProfLogLik` | `estimate_profile_loglik` | translated |
| `gridProf` | `grid_profile` | translated |
| `adapGridProf` | `adaptive_grid_profile` | translated |
| `equalTheta` | `fit_equal_theta` | translated |
| `nullTest` | `null_test` | translated |
| `rdirichlet` | `random_dirichlet` | translated |
| `simPop` | `sim_pop_sizes`, `sim_pop_equal_n` | translated |
| private `u` | `score_function` | translated/public |
| private `obsfim` | `observed_fim` | translated/public |
| private `expfim` | `expected_fim` | translated/public |
| private `loglik` | `dirmult_loglik` | translated/public |
| private `thetafim` | `theta_fim` | translated/public |
| private `mnloglik` | `multinomial_loglik` | translated/public |
| private `dbbin.ab` | internal beta-binomial PMF | translated |
| private `profU` | integrated in `estimate_profile_loglik` | translated |
| private `equalU` | integrated in `fit_equal_theta` | translated |

There is no substantive plotting code in the supplied R package. R-specific
printing, data frames, names/dimnames, lists, and S3-style presentation are not
reproduced; typed derived types and arrays are used instead.
