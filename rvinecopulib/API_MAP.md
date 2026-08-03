# API coverage map

## Direct or close numerical counterparts

| R / C++ operation | Fortran counterpart |
|---|---|
| `bicop_dist()` | `make_bicop()` |
| `dbicop()` | `bicop_model%pdf()` |
| `pbicop()` | `bicop_model%cdf()` |
| `hbicop(..., cond_var = 1)` | `bicop_model%hfunc1()` |
| `hbicop(..., cond_var = 2)` | `bicop_model%hfunc2()` |
| inverse h-functions | `hinv1()` and `hinv2()` |
| `rbicop()` | `bicop_model%simulate()` |
| `par_to_ktau()` | `bicop_model%tau()` |
| `bicop()` parameter fitting | `fit_bicop()` |
| bivariate family selection | `select_bicop()` |
| `pseudo_obs()` | `pseudo_observations()` |
| `emp_cdf()` | `empirical_cdf()` |
| `cvine_structure()` + model | `make_cvine()` / `cvine_model` |
| `dvine_structure()` + model | `make_dvine()` / `dvine_model` |
| `dvinecop()` | `cvine_model%pdf()` / `dvine_model%pdf()` |
| `rvinecop()` | `%simulate()` |
| `rosenblatt()` | `%rosenblatt()` |
| `inverse_rosenblatt()` | `%inverse_rosenblatt()` |
| Monte Carlo `pvinecop()` | `%cdf()` |
| `truncate_model()` | `%truncate()` |
| `logLik()`, AIC, BIC | `%loglik`, `%aic()`, `%bic()` |

## Adapted interfaces

- Arrays are `variables x observations`, not R's `observations x variables`.
- Pair-copula evaluation is scalar; callers loop or use array syntax at the
  application level.
- H-inverses use robust bisection for all families instead of family-specific
  C++ implementations.
- Gaussian and Student copula CDFs use deterministic Gauss-Legendre
  integration of conditional distributions.
- Family fitting uses a self-contained bounded Nelder-Mead implementation.
- Vine fitting is sequential tree-by-tree maximum likelihood.

## Not implemented

- nonparametric transformation local likelihood (`tll`)
- mixed discrete/continuous copula likelihoods and randomized discrete
  Rosenblatt transforms
- arbitrary R-vine matrices and general proximity-condition validation
- automatic Dissmann maximum-spanning-tree structure selection
- random R-vine structure generation
- threshold/truncation-level selection and mBICV
- quasi-random Halton/Sobol simulation
- parallel thread pools
- JSON/file serialization
- R marginal fitting (`kde1d`), formulas, data frames, S3 methods, summaries,
  plotting, and packaged datasets
