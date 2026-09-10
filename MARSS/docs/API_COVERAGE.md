# API coverage and compatibility notes

## Coverage basis

The upstream `NAMESPACE` exports 31 names. The computational denominator is 22 distinct exported computational functions. Nine exports are excluded under the translation rules:

- dispatch/re-export only: `accuracy`, `forecast`, `glance`, `tidy`;
- plotting/presentation: `autoplot.marssMLE`, `autoplot.marssPredict`, `CSEGriskfigure`, `CSEGtmufigure`;
- informational display: `MARSSinfo`.

Registered but unexported S3 methods and internal helpers are not counted. Mapping is **22 of 22 (100.0%)**. This is a function-mapping fraction, not a claim of complete R compatibility.

## Numerical representation

`marss_model` represents effective state-space matrices directly:

```text
x(t) = B(t) x(t-1) + U(t) + w(t),  Cov(w(t)) = Q(t)
y(t) = Z(t) x(t)   + A(t) + v(t),  Cov(v(t)) = R(t)
```

Static blocks are always present. Optional `b_t,u_t,q_t,z_t,a_t,r_t` arrays override their static blocks by time slice. `x0` and `V0` are time invariant. NaNs in `y` are missing observations.

The model validator requires consistent dimensions and symmetric positive-semidefinite covariance matrices. Filtering and smoothing use spectral pseudoinverses/pseudodeterminants where deterministic covariance directions occur.

`marss_kemcheck` adds EM-specific degeneracy checks when an affine `marss_constraints` object is supplied. It checks the spectral radius of every fully fixed `B` slice, enforces fixed `Z/A/D` rows for fixed-zero `R` diagonals, and verifies the smoothed fixed-zero-`R` fitted-value identity `E(y)=Z E(x)+A+Dd` with the upstream `1e-8` pseudoinverse tolerance. It also checks time-constant placement of fixed-zero `Q` diagonals, fixed `B` rows for zero-process-variance states, fixed `U/C` for states linked through fixed-zero `R`, potentially nonzero `Z`, and fixed-zero `Q`, fixed `U/C` rows for indirectly stochastic zero-`Q` states, time-constant potential `B` adjacency when a zero-`Q` `x0` or `U` component is estimated, and time-constant deterministic/indirectly-stochastic state classification based on the potential `B` graph. Remaining differences are chiefly exact R error-message behavior and higher-level model-construction checks performed outside the checker.

## Missing data

The filter/smoother accepts arbitrary NaN patterns. `marss_hatyt` implements the MARSS conditional-moment construction for mixed observed and missing response rows, including correlated `R`:

- conditional `E[y(t)|Y]`;
- `E[y(t)y(t)'|Y]`;
- `E[y(t)x(t)'|Y]`.

Those moments feed the `Z`, `A`, and `R` EM sufficient statistics, so `marss_kem` no longer requires complete observations.

## Time-varying matrices

Time-indexed effective `B,U,Q,Z,A,R` are used by filtering/smoothing, simulation, expected-observation moments, residual calculations, vectorization, cross-validation prediction, and innovations bootstrap calculations.

The unconstrained EM path can keep any time-varying block fixed while estimating other compatible static blocks. It deliberately rejects an attempt to estimate a block through EM when that same block is represented by a time-indexed array. Such models can instead be parameterized by affine constraints and fitted through `marss_optim_linear` or generalized-EM `marss_kem_linear`.

## Fixed/free/equality constraints

`marss_constraint_block` mirrors the MARSS vectorized model equation

```text
parameter_vector = fixed + design * beta
```

and `marss_constraints` supplies one block for `B,U,Q,Z,A,R,x0,V0,C,D,G,H,L`.

This supports:

- arbitrary fixed entries through `fixed` and zero design rows;
- independent free entries through separate design columns;
- equality/shared parameters by reusing a design column;
- fixed offsets plus free coordinates;
- designs spanning every slice of a time-indexed block.

`marss_constraint_from_affine` constructs that representation directly. `marss_constraint_from_labels` converts integer-labeled fixed/free/equality entries, and `marss_constraint_from_entries` mirrors the computational behavior of upstream `convert.model.mat()` for mixed numeric/character list matrices. Character entries may be simple shared labels or affine strings such as `2+0.5*alpha+3*beta`; repeated terms for a label are accumulated. Free-column labels are preserved, `marss_vectorized_parameter_names` produces block-prefixed names such as `B.diag`, and named start vectors can be reordered by either prefixed names or globally unambiguous raw labels.

`marss_optim_linear` maximizes the Kalman likelihood over `beta`, rejecting trial models whose covariance blocks are invalid. `marss_kem_linear` provides a generalized-EM route that holds E-step moments fixed while numerically maximizing the expected complete-data likelihood in the same beta coordinates, with an observed-likelihood nondecrease guard. Both retain the final beta vector as `fit%free_parameters`. `marss_constraints_set_start` updates the start coordinates of an existing constraint design without changing its fixed/equality structure.

`marss_hessian_linear`, `marss_fisher_i_linear`, and `marss_param_cis_linear` perform inference in those same free coordinates. Bootstrap and cross-validation refits can preserve these affine constraints and return beta-coordinate estimates.

The remaining constrained-EM difference is exact parity with the upstream package's specialized analytic M-step/update machinery and all edge-case convergence behavior, not absence of constrained EM.

## Positive-semidefinite covariance support

`Q`, `R`, `V0`, and time-indexed `Q/R` may be positive semidefinite. The native Kalman path uses a symmetric spectral pseudoinverse and log pseudodeterminant, so deterministic zero-variance state/observation directions can be propagated rather than forcing an SPD Cholesky solve. Simulation uses a symmetric eigendecomposition and therefore also supports zero eigenvalues.

The ordinary unconstrained BFGS route uses lower-Cholesky coordinates for estimated `Q`, `R`, and `V0`, so those *estimated* blocks stay positive definite. Affine constrained covariance blocks can include semidefinite structures, but invalid trial values are rejected rather than reparameterized.

## Function-level notes

The exact authoritative status, upstream source, Fortran module/source, public procedure names, and notes for each of the 22 mappings are stored in `fpm.toml` under `[[extra.translation.function]]`.

The most important differences by functional area are:

- `MARSS`, `MARSSfit`: `marss_model_spec`, `marss_build`, and `marss_from_data` translate common `form="marss"` and numerical `form="marxss"` shortcuts into the numerical model plus affine constraints and fit them directly with KEM or BFGS, including optional C/D covariates. `marss_dfa_spec`/`marss_dfa_build` provide the standard numerical DFA form and the overloaded high-level workflow. Mixed R-style numeric/character list-matrix entries and affine parameter-label expressions can be converted to the same constraints. Formula/S3 containers, factor-valued Z, arbitrary dimnames, and full R object conversion remain outside scope.
- `MARSSinits`: `marss_inits_linear` projects MARSS-style scalar/diagonal defaults into affine free coordinates, using upstream nonzero-coefficient averaging for time-varying affine blocks, and solves identified free `x0` coordinates from the first observation. `marss_inits_named` then permits partial raw or block-prefixed named beta overrides. R dimname-driven object conversion and every special list case are omitted.
- `MARSSaic`: AIC/AICc and the parametric (`AICbp`) and innovations (`AICbb`) bootstrap corrections are supplied; optional original-data `logL.star` vectors replace R object mutation. Parametric bootstrap preserves the original NaN pattern, while innovations bootstrap continues to require complete observations.
- `MARSSkem`: missing-data EM and fixed time-varying blocks are supported; `marss_kem_linear` adds generalized EM over arbitrary affine beta constraints, including designs spanning time-indexed blocks. Exact specialized upstream constrained M-step formulas and diffuse constrained EM remain different.
- `MARSSoptim`: ordinary covariance-safe BFGS plus affine f+D-beta constrained BFGS are supported; R `optim` bounds/control semantics are not reproduced.
- `MARSSkf`, `MARSSkfss`: native time-varying, missing-data, semidefinite-capable algorithms replace R/KFAS object pathways.
- `MARSSkfas`: `marss_kfas_bridge` converts the numerical model to the sibling Fortran `KFAS` representation using upstream static/time-varying shapes and `t+1` transition/disturbance indexing. It preserves raw `G/Q` disturbance coordinates, forms effective `H R H'` and `L V0 L'`, handles missing data and correlated-H LDL transformation, exposes finite and diffuse predicted/innovation covariances plus diffuse gain/rank diagnostics, reconstructs `tinitx=0` smoothed `x0/V0`, marks the first `tinitx=1` lag covariance as NaN, supports optional unstacked `return.lag.one=FALSE`, and supplies `marss_kfas_loglik` for `only.logLik`. R `SSModel`/`KFS` object-return semantics and exhaustive runtime coverage of every diffuse-rank pattern are not reproduced.
- `MARSShatyt`: conditional first/second moments used by EM and the extended supported numerical moment family are translated, including correlated-R missing conditioning; R object labeling and unsupported diffuse-object variants remain different.
- `MARSShessian`, `MARSSFisherI`, `MARSSparamCIs`: packed-parameter and constrained-beta Hessian inference is supplied. `marss_hessian_summary_linear` now provides the fitted beta mean, block-prefixed parameter names, observed-information matrix, and inverse-information covariance with `Harvey1989`, `fdHess`, and `optim` dispatch; the `optim` label shares the central finite-difference engine rather than reproducing R `optimhess` bit-for-bit. Harvey free-coordinate Fisher information is available for supported non-diffuse models, and parametric/innovations bootstrap intervals preserve affine constraints under BFGS or generalized-EM refitting. R object mutation/warning text is omitted.
- `MARSSresiduals`: basic residuals plus smoothed (`tT`-style), filtered (`tt`-style), one-step (`tt1`-style), and Harvey disturbance-smoothing routines are supplied for the supported numerical model. Exact R result-object/covariance naming conventions and diffuse Harvey handling remain different.
- `MARSSboot`, `MARSScv`: parametric/innovations bootstrap and fold/future cross-validation can preserve affine constraints, refit with BFGS or generalized EM, and return constrained beta coordinates. Bootstrap can optionally return generated data, and `marss_boot_hessian` implements `param.gen="hessian"` by drawing from the inverse observed-information covariance in raw or affine-beta coordinates. R progress/data-frame/result-object formatting and every control option are omitted.
- `MARSSinnovationsboot`: complete data remain required, matching the upstream algorithm; time-indexed model effects enter through the native filter and effective matrices.
- `MARSSsimulate`: parametric simulation supports time-varying/PSD models and optional 2-D or replicate-specific 3-D logical missingness masks corresponding to upstream `miss.loc` NaN placement; exact R RNG streams and result objects are omitted.
- `MARSSvectorizeparam`: raw active numerical blocks can be packed/unpacked, affine constrained models can be vectorized/unvectorized directly in the true free beta coordinates, custom free-column labels are retained, upstream-style block-prefixed vector names are available, and named beta inputs are reordered safely. R dimnames and arbitrary nested object fields are not reproduced.

## R interfaces intentionally outside scope

The following are intentionally not emulated as R interfaces: S3 classes/methods, general formula parsing, factor/dimname containers, plotting, printing, `data.frame`/tibble objects, progress bars, R warnings, R RNG streams, and nlme/mvtnorm object construction. The numerical list-matrix affine syntax used by `convert.model.mat()` is supported independently of R list objects. The sibling Fortran `KFAS` dependency is used numerically, not as an emulation of R KFAS objects.

The R `marxss` and `dfa` S3/list objects are not mirrored as containers. Numerically, `marss_model` exposes `C,D,G,H,L`, `marss_build` handles common `marss`/`marxss` shortcuts and covariates, `marss_dfa_spec` handles the standard DFA model vocabulary and high-level fit dispatch, and arbitrary mixed numeric/character parameter blocks can be reduced to affine constraints. Factor-valued Z, user row/column dimnames, non-affine R expressions, and automatic R object conversions remain outside scope.
