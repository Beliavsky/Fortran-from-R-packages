# MARSS - modern Fortran translation

This directory contains a modern free-form Fortran translation of the computational core of the R package **MARSS** (upstream version 3.11.10), for multivariate autoregressive state-space models.

The translation is a numerical Fortran API rather than an R runtime emulation. It preserves MARSS's state-space calculations while deliberately omitting plotting, printing, S3 dispatch, formula parsing, dimnames, and other R-specific presentation/interface behavior.

## Build

Place this directory beside the shared `rfortran-linalg` and translated `KFAS` packages in the root of `Fortran-from-R-packages`:

```text
Fortran-from-R-packages/
|-- KFAS/
|-- MARSS/
`-- rfortran-linalg/
```

Then run:

```sh
fpm build
fpm test
```

`rfortran-linalg` supplies the shared linear-algebra API and uses the repository's pinned `fortran-lapack` FPM dependency. The sibling `KFAS` translation supplies the diffuse-filter bridge. MARSS vendors neither dependency and does not link to system BLAS/LAPACK.

## Numerical model

The base model is

```text
x(t) = B(t) x(t-1) + U(t) + w(t),   w(t) ~ N(0,Q(t))
y(t) = Z(t) x(t)   + A(t) + v(t),   v(t) ~ N(0,R(t))
```

with `x0`, `V0`, observations `y`, and `tinitx=0` or `1`.

`marss_model` always contains static `B,U,Q,Z,A,R` arrays and may additionally contain `b_t,u_t,q_t,z_t,a_t,r_t`. An allocated time-indexed block overrides its static fallback. NaNs in `y` are missing observations.

The maintained numerical API supports:

- Kalman filtering, innovations likelihood, RTS smoothing, and lag-one smoothed covariances;
- static or time-indexed `B,U,Q,Z,A,R`;
- NaN-aware mixed missing/observed data;
- positive-semidefinite, including deterministic zero-variance, covariance matrices;
- Gaussian simulation with static or time-indexed model blocks and optional 2-D/3-D missingness masks matching `miss.loc` placement;
- missing-data-aware EM for unconstrained static parameter blocks;
- EM with time-varying blocks held fixed while compatible static blocks are estimated;
- covariance-safe finite-difference BFGS for ordinary numerical models;
- general affine **f + D beta** fixed/free/equality constraints through `marss_constraints`, `marss_optim_linear`, and generalized-EM `marss_kem_linear`;
- constrained-coordinate numerical Hessians, observed Fisher information, and normal confidence intervals;
- correlated-`R` conditional observation moments used by the missing-data E-step;
- fold and future cross-validation, including time-indexed prediction matrices and constraint-aware BFGS/KEM refitting;
- parametric and Stoffer-Wall innovations bootstrap paths with optional affine-constraint-preserving BFGS/KEM refitting and returned bootstrap data;
- AIC/AICc plus bootstrap AICbp/AICbb, Hessian and bootstrap parameter intervals/bias estimates, residual helpers, parameter vectorization, and basic numerical utilities.
- high-level `marss_model_spec` construction for common R `form="marss"` and numerical `form="marxss"` shortcuts, optional `C/D` covariates, and direct KEM/BFGS fitting through `marss_from_data`;
- dedicated high-level `marss_dfa_spec` construction/fitting for the standard upstream DFA model vocabulary;
- R-style mixed numeric/character list-matrix constraints through labeled entries or affine expressions, with preserved free-parameter labels;
- constraint-aware initialization and free-coordinate vectorization/unvectorization, including raw and `B.`/`Q.`-prefixed named beta ordering.

### High-level model construction

`marss_model_spec` provides common `form="marss"` and numerical `form="marxss"` structures without requiring callers to hand-build every matrix. Set `spec%form` to either `"marss"` (default) or `"marxss"`; unsupported forms such as `"dfa"` are rejected by this builder and remain available only through their dedicated numerical helpers. Supported shortcuts include `identity`, `zero`, `unconstrained`/`unequal`, `equal`, `diagonal and equal`, `diagonal and unequal`, `equalvarcov`, `onestate` for `Z`, and `scaling` for `A` when `Z` is a fixed design matrix. `G`, `H`, and `L` accept the upstream `identity`/`zero` forms.

```fortran
type(marss_model_spec) :: spec
type(marss_constraints) :: constraints
type(marss_fit_result) :: fit
integer :: info

spec%b = "diagonal and equal"
spec%a = "zero"
call marss_from_data(y, spec, 50, 1.0e-6_dp, fit, constraints, "kem", info)
```

`marss_build` exposes construction separately from fitting. Optional state and observation covariates create numerical `C*c(t)` and `D*d(t)` terms; a one-column covariate matrix is repeated over all times. `marss_from_data` accepts the same optional covariates, so supported `form="marxss"` specifications can be constructed and fitted in one call. The builder produces both the numerical `marss_model` and its affine `marss_constraints`, so the same model can flow directly into constrained EM, BFGS, bootstrap, CV, Hessian, and CI routines.

Standard dynamic-factor models use a separate `marss_dfa_spec`. It implements the computational `form="dfa"` defaults and common text-model overrides for `B`, `Q`, `A`, `R`, `x0`, `V0`, and observation-covariate `D`, including the identified lower-triangular loading pattern in `Z`. The overloaded `marss_from_data` accepts either `marss_model_spec` or `marss_dfa_spec`, so DFA models can also be constructed and fit in one call with KEM or BFGS.

`marss_inits_linear` applies MARSS-style defaults in that free-coordinate space. Scalar starts for `B`, `Q`, `R`, and `V0` are treated as diagonal targets before projection into the affine design; for time-varying affine blocks the projection now follows upstream MARSS by averaging nonzero fixed/free coefficients across time before solving. Free `x0` is solved from the first observation when the constrained system is identified. `marss_inits_named` starts from those defaults and applies a partial set of raw or block-prefixed named beta overrides, so callers can reproduce the useful computational behavior of an R `inits` list without constructing an R object.

`marss_hessian_summary_linear` mirrors the computational result of `MARSShessian`: it stores the fitted beta mean, MARSSvectorizeparam-style parameter names, observed information, and its inverse parameter covariance. `Harvey1989` is the default; `fdHess` and `optim` select the numerical information route. The `optim` label uses the same fixed-point central finite-difference engine as `fdHess`, rather than calling R's `optimhess`.

### Affine constraints

MARSS's internal constrained-model mathematics represents a parameter block as

```text
vec(parameter) = fixed + design * beta
```

The Fortran `marss_constraint_block` mirrors that representation directly. Shared columns in `design` impose equality constraints, zero rows keep fixed values, and the same representation can span all time slices of an allocated time-indexed block.

Three constructors cover progressively higher-level R model-matrix conventions:

- `marss_constraint_from_affine` accepts the already reduced `fixed`, `design`, and beta-start arrays directly;
- `marss_constraint_from_labels` uses integer labels, where zero means fixed and repeated positive labels impose equality constraints;
- `marss_constraint_from_entries` mirrors the computational part of upstream `convert.model.mat()` for mixed numeric/character list matrices. Character entries can be simple shared labels or affine expressions such as `2+0.5*alpha+3*beta`; repeated terms are accumulated and first-occurrence parameter names are preserved.

`marss_free_parameter_names` returns raw free-column labels. `marss_vectorized_parameter_names` adds the upstream `MARSSvectorizeparam()` block prefixes such as `B.` and `Q.`, while `marss_reorder_free_parameters` and `marss_constraints_set_start_named` accept those prefixed names or globally unambiguous raw labels.

`marss_optim_linear` fits these coordinates directly by BFGS, while `marss_kem_linear` performs a generalized-EM fit by maximizing the expected complete-data log likelihood in the same beta coordinates. Both store the final `beta` vector in `fit%free_parameters`. `marss_constraints_set_start` can copy a fitted beta vector back into an existing constraint design without changing fixed/equality structure. The companion `marss_hessian_linear`, `marss_fisher_i_linear`, and `marss_param_cis_linear` operate in that same free-coordinate space.

The generalized-EM route is numerical rather than a verbatim port of every specialized analytic constrained M-step in R MARSS, so iteration paths and edge-case convergence can differ even when the fitted constrained likelihood problem is equivalent.

## Minimal filtering example

```fortran
program demo
   use marss_api
   implicit none

   type(marss_model) :: model
   type(marss_kf_result) :: kf

   allocate(model%y(1, 4), model%b(1, 1), model%u(1), model%q(1, 1))
   allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
   model%y(1, :) = [0.2_dp, 0.4_dp, 0.3_dp, 0.7_dp]
   model%b = 0.85_dp
   model%u = 0.05_dp
   model%q = 0.08_dp
   model%z = 1.0_dp
   model%a = 0.0_dp
   model%r = 0.12_dp
   model%x0 = 0.0_dp
   model%v0 = 0.5_dp
   model%tinitx = 0

   call marss_kfss(model, kf)
   if (.not. kf%ok) error stop "Kalman filter failed"
end program demo
```

See `example/kalman_example.f90`, `example/em_example.f90`, `example/constraints_example.f90`,
`example/high_level_model.f90`, `example/dfa_high_level.f90`, and `example/named_constraints.f90`.

## Remaining compatibility differences

The largest remaining differences are now concentrated in areas that are either R-specific or specialized algorithms:

- R formula objects, S3 classes/methods, dimnames, factors, warnings, and result-object formatting are not reproduced.
- Standard numerical DFA construction/dispatch is now available through `marss_dfa_spec`, and mixed list-matrix labels/affine expressions can be converted into `f + D beta` blocks. Remaining construction differences are the R list/factor container itself, factor-valued `Z`, arbitrary expression parsing beyond upstream's affine `a+b*p` convention, user row/column dimnames, and automatic conversion between R `marxss`/`marss` objects.
- Constrained generalized EM is implemented in affine beta coordinates, but it does not reproduce every specialized analytic upstream constrained M-step or all degeneracy/convergence diagnostics exactly.
- Diffuse models dispatch through the sibling `KFAS` bridge. The bridge now mirrors upstream static-versus-time-varying storage, MARSS/KFAS transition indexing, raw `G/Q` disturbance coordinates, effective `H R H'` and `L V0 L'` covariances, missing data, exact finite/diffuse covariance decomposition, `tinitx=0/1` initial smoothing, `return.lag.one`, and `only.logLik`. Remaining differences are R `SSModel`/`KFS` object-return semantics and runtime validation over every possible diffuse-rank pattern.
- Harvey free-coordinate Fisher information and `tT`/`tt`/`tt1`-style plus Harvey residual routines are implemented for the supported non-diffuse numerical model, but exact R object/covariance conventions and diffuse variants are not complete.
- Bootstrap and cross-validation preserve affine constraints and can refit with BFGS or generalized EM; Hessian-normal parameter generation is also available. R RNG streams, progress/result-object formatting, and all R optimizer-control semantics are not reproduced.
- Free-parameter vectorization now preserves custom free-column labels, produces upstream-style block-prefixed vector names, and supports name-based reordering. R dimnames and arbitrary nested list/object fields remain outside the numerical API.

`marss_kemcheck` now implements the principal fixed-zero variance and graph restrictions used by KEM: fully fixed `B` slices must have spectral radius at most one, zero observation-error rows require fixed `Z/A/D`, and fixed-zero `R` rows must satisfy the smoothed fitted-value identity `E(y)=Z E(x)+A+Dd` within the upstream `1e-8` pseudoinverse tolerance. It also requires time-invariant zero-process-variance placement, fixed `B` rows for zero-process-variance states, fixed `U/C` for states linked through `R=0`, potentially nonzero `Z`, and `Q=0`, fixed `U/C` on indirectly stochastic zero-`Q` states, time-invariant potential `B` adjacency when deterministic-state `x0` or `U` is estimated, and time-invariant deterministic versus indirectly stochastic state classification. Remaining KEM differences are chiefly exact R diagnostic-message behavior and model-construction checks that upstream performs outside `MARSSkemcheck`.

These limitations are compatibility differences, not additional untranslated exported-function counts.

## Translation coverage

Package status: **substantial**.

**22 of 22 (100.0%)** exported computational R functions have a meaningful Fortran mapping. This percentage measures **mapped computational functions**, not complete R interface or behavioral compatibility. Plotting, printing, formatting, generic re-exports/dispatch-only functions, informational displays, datasets, and presentation-only exports are excluded from the denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `MARSS` | `marss_em`, `marss_builder`, `marss_dfa_mod`, `marss_workflow` | `marss`, `marss_build`, `marss_dfa_build`, `marss_from_data` | substantial |
| `MARSSaic` | `marss_analysis`, `marss_bootstrap` | `marss_aic`, `marss_bootstrap_aic` | substantial |
| `MARSSboot` | `marss_bootstrap` | `marss_boot`, `marss_boot_hessian` | substantial |
| `MARSScv` | `marss_cv_mod` | `marss_cv` | substantial |
| `MARSShessian` | `marss_analysis`, `marss_constraints_mod`, `marss_hessian_summary_mod` | `marss_hessian`, `marss_hessian_linear`, `marss_hessian_summary_linear` | substantial |
| `MARSSinits` | `marss_model_ops`, `marss_initialization` | `marss_inits`, `marss_inits_linear`, `marss_inits_named` | substantial |
| `MARSSkem` | `marss_em`, `marss_constrained_em` | `marss_kem`, `marss_kem_linear` | substantial |
| `MARSSkemcheck` | `marss_kemcheck_mod` | `marss_kemcheck` | partial |
| `MARSSkf` | `marss_kalman` | `marss_kf` | substantial |
| `MARSShatyt` | `marss_analysis` | `marss_hatyt` | substantial |
| `MARSSkfss` | `marss_kalman` | `marss_kfss` | substantial |
| `MARSSkfas` | `marss_kfas_bridge` | `marss_kfas`, `marss_kfas_loglik` | substantial |
| `MARSSoptim` | `marss_optim_mod`, `marss_constraints_mod` | `marss_optim`, `marss_optim_linear` | substantial |
| `MARSSparamCIs` | `marss_analysis`, `marss_constraints_mod`, `marss_bootstrap` | `marss_param_cis`, `marss_param_cis_linear`, `marss_bootstrap_param_cis` | substantial |
| `MARSSresiduals` | `marss_analysis`, `marss_residuals_full` | `marss_residuals`, `marss_residuals_smoothed`, `marss_residuals_filtered`, `marss_residuals_one_step`, `marss_residuals_harvey` | substantial |
| `MARSSsimulate` | `marss_simulation_mod` | `marss_simulate` | substantial |
| `MARSSFisherI` | `marss_analysis`, `marss_constraints_mod`, `marss_harvey_info` | `marss_fisher_i`, `marss_fisher_i_linear`, `marss_fisher_i_harvey_linear` | substantial |
| `MARSSvectorizeparam` | `marss_analysis`, `marss_constraints_mod` | `marss_vectorizeparam`, `marss_unvectorizeparam`, `marss_vectorize_free`, `marss_unvectorize_free`, `marss_free_parameter_names`, `marss_vectorized_parameter_names`, `marss_reorder_free_parameters` | substantial |
| `zscore` | `marss_utils` | `zscore` | complete |
| `ldiag` | `marss_utils` | `ldiag` | substantial |
| `MARSSinnovationsboot` | `marss_innovations` | `marss_innovations_boot` | substantial |
| `MARSSfit` | `marss_em`, `marss_dfa_mod`, `marss_workflow` | `marss_fit`, `marss_dfa_fit_spec`, `marss_from_data` | substantial |

Untranslated computational exports: none.

The authoritative machine-readable mapping and per-function notes are in `fpm.toml`. `docs/API_COVERAGE.md` gives additional compatibility detail.

## Tests

`test/test_marss.f90` is deterministic. It includes KFAS conversion/runtime parity cases for static and time-varying systems, raw `G/Q` and `H/R`/`L/V0` loading coordinates, correlated observation covariance requiring LDL transformation, fully missing time points, rank-deficient process noise, singular observation noise, diffuse-rank diagnostics, `tinitx=0` initial smoothing, stacked versus unstacked lag-one models, and log-likelihood-only evaluation. It also exercises high-level MARSS/MARXSS and DFA construction, mixed list-matrix label/expression conversion, raw/prefixed free-parameter naming and reordering, affine default initialization, raw and free-coordinate vectorization, direct KEM/BFGS workflow fitting, static and time-varying filtering, fully deterministic PSD cases, mixed missing observations, correlated-`R` missing-data moments, missing-data EM, fixed-time-varying EM, simulation, vectorization, unconstrained and affine-constrained BFGS, constrained Hessian/CIs, both explicit and default cross-validation folds including affine-constrained BFGS/KEM refits, constraint-preserving parametric/innovations bootstrap with returned beta coordinates/data, Hessian-normal parameter generation, bootstrap AIC and bootstrap parameter-CI/bias routines, and the earlier fixed numerical Kalman reference, and KEM fixed-zero degeneracy diagnostics, including fixed-`B` unit-circle checks.

## License and provenance

The upstream package declares `GPL-2`. This derivative translation is distributed under GPL-2; see `LICENSE` and `NOTICE.md`. Upstream metadata, citation information, NAMESPACE, and R source retained for provenance are under `upstream/`. No source from `rfortran-linalg`, KFAS, mvtnorm, nlme, BLAS, LAPACK, ARPACK, or other dependency packages is copied into MARSS.
