# R mgcv to Fortran API map

## Direct or close numerical translations

| R mgcv concept/export | Fortran API | Notes |
|---|---|---|
| `gam.fit`, core `gam` fitting | `gam_fit` | Dense penalized IRLS; caller supplies the design and penalties. |
| `magic` | `magic_fit` | Dense Gaussian generalized-ridge fit. |
| GCV/UBRE smoothing selection | `method_gcv`, `method_ubre` | Deterministic log-smoothing-parameter coordinate search. |
| REML smoothing selection | `method_reml` | Dense REML-like objective; not the full derivative engine in R mgcv. |
| `Predict.matrix*` | `predict_smooth`, `predict_gam` | Typed smooth metadata replaces S3 dispatch. |
| cubic regression smooth | `construct_cr_smooth` | Natural cubic basis with integrated squared second-derivative penalty. |
| `ps` smooth | `construct_ps_smooth` | B-spline basis plus finite-difference penalty. |
| cyclic smooth | `construct_cyclic_smooth` | Periodic Fourier basis with fourth-order frequency penalty. |
| thin-plate smooth | `construct_tp_smooth_1d`, `construct_tp_smooth_2d` | Low-rank radial approximation; see porting notes. |
| random-effect smooth | `construct_random_effect` | One-hot design with identity penalty. |
| `tensor.prod.model.matrix` | `tensor_product_model_matrix`, `tensor_smooth` | Two-margin row tensor product. |
| `tensor.prod.penalties` | `tensor_product_penalties` | Kronecker-sum margin penalties. |
| `place.knots` | `place_knots` | Type-7 quantile placement. |
| `pcls` | `pcls_fit` | Dense projected penalized least squares. |
| `mono.con` | `monotonicity_constraints` | First-difference constraints on ordered fitted values/coefficient grids. |
| `bandchol`, `trichol` | `band_cholesky`, `tri_cholesky` | Native band-storage/tridiagonal algorithms. |
| `Rrank` | `rrank` | Symmetric-product eigenvalue rank estimate. |
| `mroot`, `mini.roots` | `mroot`, `mini_root` | Symmetric eigen roots. |
| `slanczos` | `slanczos` | Dense eigensystem replacement rather than iterative Lanczos. |
| `notExp`, `notLog` | `not_exp`, `not_log` | Same piecewise transformations. |
| `notExp2`, `notLog2` | `not_exp2`, `not_log2` | Same bounded-log-range transformation. |
| `dpnorm` | `dpnorm` | Cancellation-aware normal interval probability. |
| `rmvn`, `dmvn` | `rmvn`, `dmvn_log` | Native dense covariance routines. |
| `r.mvt`, `d.mvt` | `rmvt`, `dmvt_log` | Native multivariate t routines. |
| `rTweedie` | `rtweedie` | Compound Poisson-gamma case, `1 < p < 2`. |
| `psum.chisq` | `weighted_chisq_cdf` | Numerical characteristic-function inversion. |
| `XWXd`, `XWyd` | `xwxd`, `xwyd` | Dense model-matrix equivalents. |
| `Xbd` | `xbd` | Dense matrix-vector product. |
| `diagXVXd`, `ijXVXd` | `diag_xvxd`, `ij_xvxd` | Dense covariance contractions. |
| `uniquecombs` | `unique_rows` | Numeric matrices with optional tolerance. |
| `exclude.too.far` | `exclude_too_far` | Scaled Euclidean nearest-data test. |
| `null.space.dimension` | `null_space_dimension` | Polynomial null-space dimension. |
| `gamSim` | `gam_sim` | Examples 1 and 2; normal, Poisson, and binary. |

## Families implemented

`family_t` supports:

- `family_gaussian`
- `family_binomial`
- `family_poisson`
- `family_gamma`
- `family_inverse_gaussian`
- `family_negative_binomial` with fixed `theta`
- `family_tweedie` with fixed `p`

## Deliberately omitted or not yet translated

The following are substantial R subsystems rather than standalone numerical
routines and are not represented as if they were complete:

- Formula parsing, model frames, contrasts, factors, S3 methods, printing,
  summaries, data-frame/`xts` metadata, and missing-data class restoration.
- Plotting (`plot.gam`, `vis.gam`, `qq.gam`, polygon/soap plotting).
- `bam`'s compressed/discrete sparse engine, sparse QR/Cholesky paths,
  OpenMP scheduling, and AR residual machinery. Dense cross-product
  equivalents are supplied instead.
- Full `magic`/`gdi` derivative calculations, NCV, Fellner-Schall updates,
  exact outer Newton optimization, and exact mgcv REML/ML constants.
- `gamm`/`nlme`, `pdTens`, `pdIdnot`, correlation structures, and mixed-model
  object conversion.
- Adaptive, Duchon, soap-film, spherical, MRF, Gaussian-process, SCAD,
  shrinkage, factor-smooth, and shape-constrained smooth classes.
- General/extended family systems (`gamlss`, `gfam`, `gevlss`, `gaulss`,
  `twlss`, `ziplss`, `multinom`, `mvn`, `cox.ph`, ordinal models, etc.).
- JAGS, INLA, full Bayesian Laplace approximation, posterior GAM MCMC,
  survival-specific baseline-hazard utilities, and model diagnostics.

The original sources for all of these remain in `original/` for provenance and
future incremental translation.
