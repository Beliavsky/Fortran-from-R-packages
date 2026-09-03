# API coverage

This document records computational parity against MCMCglmm 2.36. "Implemented"
means there is a maintained Fortran numerical API with deterministic regression
coverage. "Partial" means the numerical algorithm is translated but R formula,
S3/coda, R sparse-object, or data-frame orchestration is intentionally outside
the Fortran surface. Native Fortran CSR design inputs are supported separately.

## Main MCMC engines

| Upstream area | Status | Maintained Fortran API / notes |
| --- | --- | --- |
| Gaussian multivariate mixed models | Implemented dense core | `gaussian_mixed_mcmc` jointly samples fixed/random effects and `G`/`R`. |
| Multiple random-effect covariance structures | Implemented | `multi_term_gaussian_mixed_mcmc` and the multi-term non-Gaussian/grouped engines assign concatenated random-effect columns to independent `G` blocks. |
| Heterogeneous scalar responses | Implemented dense core | `heterogeneous_family_mixed_mcmc` and `heterogeneous_multi_term_mixed_mcmc` mix Gaussian and scalar non-Gaussian traits in a full latent residual covariance block, with response masks and imputation. |
| Grouped multinomial families | Implemented dense core | `multinomial_family_mixed_mcmc` and `multinomial_multi_term_mixed_mcmc` cover multinomial, `ztmb`, and zero-truncated multinomial blocks. |
| Zero-inflated/hurdle/zero-altered families | Implemented dense core | `two_part_mixed_mcmc` and `two_part_multi_term_mixed_mcmc` use correlated two-process latent models. |
| Ordered-probit family 14 | Implemented | `ordinal_native_mixed_mcmc`; native observation layer, liability updates, adaptive cutpoint updates, and the upstream binary slice optimization when requested. |
| Threshold family 20 | Implemented | `threshold_cutpoint_mixed_mcmc`; truncated-normal liabilities and MCMCglmm-style sequential adaptive cutpoints. Fixed-cutpoint `ordinal_probit_mixed_mcmc` remains available. |
| Missing responses | Implemented across maintained engines | Optional masks omit the observation likelihood and impute through the latent model for heterogeneous, grouped, two-process, ordinal, and threshold paths. |
| Native binary slice updates | Implemented | `binary_slice_liability_update` reproduces the three upstream univariate binary slice transformations used for families 3, 14, and 22. |
| Adaptive liability proposals | Implemented numerical core | `optimal_acceptance_ratio`, `adaptive_mh_observe`, `adaptive_mh_decay`, and `adaptive_mh_finalize` reproduce the native target/decay recursions; maintained family engines can use them where applicable. |
| `theta_scale` | Implemented Gaussian multi-term engine | `theta_scale_gaussian_mixed_mcmc` plus `theta_scale_multivariate_conditional`; numeric design masks replace R formula bookkeeping. |
| Parameter expansion | Implemented in several dense engines | Single-`G` Gaussian, multi-`G` Gaussian, and heterogeneous scalar multi-`G` latent parameter-expanded chains are implemented. Standalone `parameter_expansion_conditional`, `apply_parameter_expansion`, and `expanded_covariance` are public. Grouped/two-process/ordinal parameter expansion is not yet unified. |
| Path/SIR structural coefficients | Implemented Gaussian multi-term engine | `structural_gaussian_multi_term_mcmc` applies `Lambda = I - sum(lambda_k B_k)`, samples transformed mixed models, and updates structural coefficients using the Jacobian-aware MH target. `sir_matrix` supplies the numerical SIR interaction matrix. |
| Measurement error `me(type="dberkson")` | Implemented in scalar chain | `categorical_measurement_error_update` samples latent categories from class priors times MVN residual likelihoods and is integrated into the heterogeneous scalar route of `mcmcglmm_fit_numeric`. R formula/group construction is skipped; other `me` types are also unimplemented upstream. |
| `covu` joint random/residual covariance | Implemented specialized Gaussian dense engine | `covu_gaussian_mixed_mcmc` samples a joint G-R covariance and reports marginal `G`, Schur-complement conditional `R`, and residual-on-random regression. `joint_gr_decompose`, `joint_gr_compose`, and `joint_gr_covariance_update` expose the matrix core. Upstream itself disallows `ginverse` and prediction/simulation for `covu`; the Fortran API does not invent those unsupported combinations. |
| Sparse production solver | Partial implementation | `mcmcglmm_sparse_matrix` provides canonical one-based CSR algebra. The heterogeneous scalar route forms predictors and the stacked design crossproduct from paired CSR `X`/`Z`, assembles its joint coefficient precision as CSR coordinate sums, and uses independent sparse Cholesky with deterministic reverse Cuthill-McKee ordering. A column-assisted row accumulator avoids storing all duplicate observation-level crossproduct contributions. Ordering and symbolic analysis are cached and safely rebuilt after structural expansion. Other routes currently materialize CSR designs. Modified upstream CSparse is deliberately not copied or vendored. |
| Common numeric orchestration | Implemented foundation | `mcmcglmm_fit_numeric` validates and routes `mcmcglmm_numeric_model`, `mcmcglmm_numeric_prior`, and `mcmcglmm_control` through nine engine classes: heterogeneous scalar, parameter-expanded scalar, two-process, grouped multinomial, ordinal, threshold, `theta_scale`, structural path/SIR, and `covu`. It accepts unambiguous dense or CSR designs. Covariance update masks, modes, splits, fixed blocks, and antedependence options use one normalized configuration where supported by the underlying engine. |
| Single all-features orchestration equivalent to R `MCMCglmm()` | Partial | The common numeric entry point now reaches the principal translated engines, but it does not combine unrelated response-block classes in one joint sparse chain. Discrete Berkson measurement error is integrated with the scalar route; arbitrary cross-feature combinations and formula-driven R object construction remain interface work. |

## Native family-code coverage

MCMCglmm's native C++ switch uses integer family codes. The maintained Fortran
likelihood and response-generation layers cover every computationally active
branch; branches that upstream leaves empty or explicitly errors are not invented.

| Code | Upstream branch | Status |
| ---: | --- | --- |
| 1 | Gaussian | Implemented |
| 2 | Poisson | Implemented |
| 3 | Multinomial | Implemented grouped; binary/binomial scalar representation and native binary slice update also available |
| 4 | shape-one Weibull kernel (`notyet_weibull` label) | Implemented native shape-one kernel |
| 5 | Exponential | Implemented |
| 6 | Censored Gaussian | Implemented |
| 7 | Censored Poisson | Implemented |
| 8 | shape-one censored Weibull kernel (`notyet_cenweibull` label) | Implemented native shape-one censored kernel |
| 9 | Censored Exponential | Implemented |
| 10 | Zero-inflated Gaussian | Upstream native switch is empty; not invented |
| 11 | Zero-inflated Poisson | Implemented |
| 12 | Zero-inflated Weibull | Upstream native switch is empty; not invented |
| 13 | Zero-inflated Exponential | Upstream native switch is empty; not invented |
| 14 | Ordered multinomial probit | Implemented, including adaptive cutpoints and binary slice path |
| 15 | Hurdle Poisson | Implemented |
| 16 | Zero-truncated Poisson | Implemented |
| 17 | Geometric, support starting at zero | Implemented with upstream link sign convention |
| 18 | Zero-altered Poisson | Implemented |
| 19 | Zero-inflated Binomial | Implemented |
| 20 | Threshold | Implemented, including adaptive cutpoints |
| 21 | `zitobit` | Upstream explicitly raises "not yet implemented"; not invented |
| 22 | Nonzero Binomial | Implemented, including native binary slice path |
| 23 | Noncentral scaled-t | Implemented; deterministic density regressions agree with SciPy reference values |
| 24 | Mean-shifted scaled-t | Implemented |
| 25 | Hurdle Binomial | Implemented |
| 26 | Zero-truncated multiple Bernoulli (`ztmb`) | Implemented |
| 27 | Zero-truncated Multinomial | Implemented |

## Covariance update-code coverage

The numerical dispatcher `covariance_update_dispatch` implements the covariance
operations behind native update codes 0 through 6. Modes 0-6 are integrated into
the multi-term Gaussian, heterogeneous scalar, two-process, and multinomial
multi-`G` engines. Mode 5 additionally needs the level-by-trait locations used by
the antedependence likelihood, which the high-level engines provide from random
effects or residuals.

| Code | Meaning | Fortran status |
| ---: | --- | --- |
| 0 | Fixed covariance | Integrated; supplied covariance remains unchanged |
| 1 | Unstructured inverse-Wishart | `inverse_wishart_sample`, `riw_mcmcglmm`, dispatcher |
| 2 | Conditioned/block-constrained IW | `conditioned_covariance_update`, conditioned `rIW`, dispatcher |
| 3 | Correlation/fixed-variance correlation | `correlation_structure_update`, dispatcher |
| 4 | Correlation-constrained submatrix | `correlation_submatrix_update`, dispatcher |
| 5 | Antedependence | `ante_covariance_samples`, dispatcher and high-level sampler routes |
| 6 | Unstructured plus identity direct sum | `identity_direct_sum_update`, dispatcher |

`covu` uses its own joint G-R routing because its sufficient statistics and Schur
structure differ from ordinary independent G/R updates.

## Exported numerical APIs

| R API / area | Status | Fortran equivalent / notes |
| --- | --- | --- |
| `MCMCglmm` | Substantial computational core | Focused and unified samplers described above; numeric orchestration accepts dense or CSR designs, while formula/S3/coda construction and a sparse factorization engine remain outside the Fortran API. |
| `inverseA` pedigree | Implemented dense core | `pedigree_relationship`, `pedigree_inverse` |
| `inverseA` phylogeny | Implemented dense core | `phylogenetic_precision`, reusing sibling `ape` |
| `rbv` pedigree | Implemented | `breeding_values_pedigree` |
| `rbv` phylogeny | Implemented numerical core | `breeding_values_phylo`, reusing sibling `ape` |
| `prunePed` | Implemented ancestry core | `prune_pedigree_mask`; data-frame reconstruction and `make.base` interface omitted |
| `rIW` | Implemented | `riw_mcmcglmm`, `riw_mcmcglmm_conditioned` |
| `rtnorm` | Implemented | `truncated_normal_sample`, including stable finite-tail sampling |
| `rtcmvnorm` | Implemented | `truncated_conditional_mvn_sample` |
| `dcmvnorm` | Implemented | `conditional_mvn_parameters`, `conditional_mvn_log_density` |
| `posterior.cor` | Implemented | `posterior_correlations` |
| `posterior.inverse` | Implemented | `posterior_inverses` |
| `posterior.evals` | Implemented | `posterior_eigenvalues` |
| `posterior.ante` | Implemented | `ante_parameters` |
| `posterior.mode` | Implemented numerical KDE | `posterior_modes`; Gaussian KDE using the R `bw.nrd0` bandwidth rule and default `adjust=0.1` concept |
| `ante`/`rante` covariance core | Implemented | `ante_covariance_samples` |
| `commutation` | Implemented | `commutation_matrix` |
| `Tri2M` | Implemented | `triangle_to_matrix`, `matrix_to_triangle` |
| `Ptensor` | Implemented | `central_moment_tensor` |
| `kunif` | Implemented | `uniform_central_moment` |
| `KPPM` | Implemented | `symmetrizer_matrix` |
| `knorm` | Implemented | `normal_moment_matrix` |
| `krzanowski.test` | Implemented numerical core | `krzanowski_compare` |
| `list2bdiag` | Implemented numerical core | `block_diagonal` |
| `gelman.prior` | Implemented numerical core | `gelman_prior_design`; formula/model-matrix construction omitted |
| `Ddivergence` | Implemented numerical core | `d_divergence_mc` |
| `mult.memb` | Implemented numerical core | `multiple_membership_design` |
| `spl` | Implemented LRTP core | `spline_lrtp`, including explicit knots and R type-7 default quantile knots |
| `path` | Implemented numerical core and Gaussian chain | `path_matrix`, structural transformation helpers, `structural_gaussian_multi_term_mcmc` |
| `sir` | Implemented numerical matrix core | `sir_matrix`; formula/model-matrix construction is interface code |
| `me(type="dberkson")` | Implemented in scalar chain | `categorical_measurement_error_update` and heterogeneous scalar orchestration; R formula/group attributes omitted |
| `buildV` | Implemented dense multi-term numerical core | `multi_term_build_v` assembles `R kron I + sum G_t kron (Z_t A_t Z_t')` and marginal variances |
| `predict.MCMCglmm` | Partial computational core | `posterior_linear_predictor` supports marginalizing selected random terms; `scalar_response_expectation` integrates latent Gaussian variance for supported scalar links. `newdata`, factor matching, interval/coda/S3 handling remain R interface work. |
| `simulate.MCMCglmm` | Substantial numerical core | Family response generators and `simulate_multi_term_gaussian_latent` draw new multi-term random effects/residuals from numeric designs and posterior covariance draws. R `newdata`/formula matching remains interface code. |
| `pkk` | Implemented | `pkk_probability` |
| `Dtensor`, `evalDtensor`, `Dexpressions` | Skipped | Symbolic R expression differentiation, not a numerical Fortran kernel |
| `sm2asreml` | Skipped | External-package/R object conversion helper |
| `at.set`, `at.level` | Skipped | R formula-construction helpers |
| `summary.MCMCglmm` | Partial via lower-level statistics | Many covariance/posterior transforms are present; coda/S3 table formatting is not translated |
| `residuals.MCMCglmm` | Not invented | Upstream method itself reports that residuals are not yet implemented |
| Plot methods | Skipped | Plotting/interface code |

## Additional maintained numerical pieces

- explicit deterministic `rng_state` with state-passed uniform, normal, gamma,
  chi-square, exponential, and Poisson variate generation;
- inverse-Wishart and conditioned inverse-Wishart samplers;
- MVN density, covariance/precision samplers, and conditional/truncated MVN tools;
- covariance-to-correlation, Kronecker, commutation, block-diagonal, and symmetric
  eigenvalue helpers;
- standalone and integrated parameter-expansion conditionals and output rescaling;
- `theta_scale` conjugate Gaussian conditional and multi-term Gaussian chain;
- exact MCMCglmm acceptance-target/adaptive-MH recursions;
- native binary slice updates used by upstream families 3, 14, and 22;
- joint G-R Schur decomposition/composition/update machinery and a specialized
  Gaussian `covu` chain;
- structural path/SIR matrices, Jacobian-aware structural likelihood, and an
  integrated Gaussian path/SIR MCMC engine;
- discrete Berkson measurement-error category updating;
- canonical one-based CSR storage, construction, validation, dense conversion,
  products, transpose products, sparse crossproducts, and sparse-transpose-times-
  sparse products;
- a CSR-native heterogeneous scalar path for predictor evaluation and coefficient
  normal-equation formation, including discrete Berkson measurement error;
- independent natural-order sparse Cholesky factorization, triangular precision
  solves, and Gaussian precision sampling;
- direct COO/CSR assembly of the scalar route's joint coefficient precision,
  with dense-versus-sparse conditional-mean regression coverage;
- CSR stacked-design crossproducts with bounded duplicate-contribution workspace,
  and deterministic reverse Cuthill-McKee symmetric ordering;
- separate symbolic/numerical sparse Cholesky phases with per-chain cache reuse
  and structural-expansion invalidation;
- posterior predictive latent Gaussian simulation plus response-generation
  transformations for the maintained family codes;
- dense multi-term marginal covariance assembly and response-scale posterior
  expectations for the main scalar links.
- deterministic orchestration tests for all nine routed engine classes, including
  draw-for-draw agreement between routed and direct scalar, ordinal, threshold,
  `theta_scale`, structural path/SIR, and `covu` entry points, plus exact dense/CSR
  input parity for the scalar route.

## Important differences and remaining targets

1. The typed orchestrator accepts dense or canonical CSR design inputs. The
   heterogeneous scalar route retains paired CSR designs for predictor and
   normal-equation operations and uses CSR posterior-precision assembly plus
   sparse Cholesky. Its stacked crossproduct is CSR and the factorization uses
   reverse Cuthill-McKee ordering with cached symbolic analysis. Its
   column-assisted accumulator does not retain duplicate observation-level pair
   contributions, but other routes still materialize CSR inputs. The package does
   not copy modified CSparse. Very large models still need production-scale
   benchmarking, stronger ordering strategies, and broader sparse latent-state
   integration for comparable scaling.
2. Randomness is supplied through an explicit `rng_state`, not R's global RNG.
   Identical Fortran seeds are reproducible within this package, but draws are not
   expected to be bit-identical to R/MCMCglmm.
3. APIs accept numeric design matrices, covariance structures, pedigree indices,
   and translated `ape` tree objects. R formula parsing, factors/contrasts,
   environment capture, S3/coda objects, and plotting are intentionally omitted.
4. `mcmcglmm_fit_numeric` now provides a common typed entry point for nine major
   engine classes. The largest remaining computational target is a shared sparse
   latent-state architecture that can combine arbitrary scalar, grouped, ordinal,
   threshold, and two-process blocks with special features in one joint chain.
5. Parameter expansion is not yet integrated into every grouped/ordinal/two-process
   engine, structural/path updates are currently Gaussian, and the specialized
   `covu` engine is Gaussian only, matching the tractable native numerical core
   without inventing unsupported upstream prediction/simulation behavior.
