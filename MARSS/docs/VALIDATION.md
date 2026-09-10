# Validation notes

The maintained deterministic suite is `test/test_marss.f90`. Runnable examples are under `example/`, including `high_level_model.f90` for specification-driven construction/fitting, `dfa_high_level.f90` for DFA dispatch, and `named_constraints.f90` for R-style affine list-matrix labels.

The intended repository verification commands are:

```sh
fpm build
fpm test
fpm clean --all
```

The generation environment had GNU Fortran 14.2 but no FPM executable. Attempts to obtain FPM were blocked by the environment's outbound download restrictions, so the literal FPM commands could not be executed during generation.

Instead, every maintained MARSS source file, the tests, and all examples were compiled directly with GNU Fortran using strict checking. The final validation build uses:

```text
-std=f2018 -O2 -Wall -Wextra -Wpedantic -fcheck=all
```

For direct compiler validation only, external temporary stand-ins for the exact `rfortran-linalg` interfaces and the compile-time `KFAS` interfaces imported by MARSS are kept outside this package tree. They are not included in the package. The translation's FPM manifest declares the sibling dependencies:

```toml
rfortran-linalg = { path = "../rfortran-linalg" }
KFAS = { path = "../KFAS" }
```

No system BLAS/LAPACK link directives or copied dependency implementations are present.

Because the temporary KFAS stand-in is interface-only, the direct validation copy skips the runtime KFAS calls while retaining all KFAS model-conversion tests. The official packaged `test/test_marss.f90` keeps those runtime tests unchanged. They compare the sibling backend against the native smoother for a multivariate time-varying model with correlated observation covariance, missing data, rank-deficient process noise, singular observation noise, and `tinitx=0`; they also exercise diffuse rank/finite-vs-infinite covariance diagnostics, `return.lag.one=FALSE`, and `only.logLik`. A repository FPM run with the actual sibling `KFAS` package remains required for literal backend runtime verification in this environment.

The deterministic suite covers, among other cases:

- a fixed numerical Kalman/smoother reference with one missing observation;
- KFAS conversion parity for static/time-varying shapes, MARSS-to-KFAS transition indexing, raw `G/Q`, effective `H R H'`, `L V0 L'`, stacked/unstacked lag-one models, missingness, and `tinitx=0/1` initialization;
- official sibling-KFAS runtime parity tests for correlated-H LDL transforms, a fully missing time point, semidefinite process noise, singular observation covariance, diffuse diagnostics, lag-one availability, and log-likelihood-only evaluation;
- time-varying `B/U/Q/Z/A/R` in a fully deterministic PSD model;
- static and time-indexed parameter vectorization round trips;
- common `marss_model_spec` `marss`/`marxss` shortcut construction, including `onestate`/`scaling`, numerical C/D covariates, and rejection of unsupported builder forms;
- high-level `marss_from_data` fitting through both generalized EM and BFGS;
- affine free-coordinate vectorization/unvectorization and projection from numerical model values;
- constraint-aware MARSS-style initialization, including diagonal scalar targets and identified free-`x0` least-squares starts;
- correlated-`R` conditional moments for a mixed observed/missing response vector;
- missing-data EM, including recovery of scalar `A` and `R` observed-data MLEs;
- EM with fixed time-varying blocks and an estimated static block;
- deterministic Gaussian simulation;
- ordinary covariance-safe BFGS likelihood improvement;
- affine f+D-beta equality-constrained BFGS;
- constrained beta-coordinate Hessian and confidence intervals;
- explicit-fold and default-fold cross-validation, including affine-constrained BFGS and generalized-EM refits with returned fold beta coordinates;
- AIC/AICc, parametric/innovations bootstrap AIC corrections, residuals, Hessian/Fisher information, normal confidence intervals, and parametric/innovations bootstrap confidence intervals with bias estimates;
- parametric and innovations-bootstrap reproducibility, direct bootstrap-data return, affine-constraint-preserving BFGS/KEM refits, and constrained beta-coordinate outputs;
- KEM fixed-zero R/Q degeneracy restrictions, including fixed-`B` spectral-radius checks (real and complex eigenvalue cases), fixed-zero-`R` smoothed fitted-value consistency, time-varying Q zero-placement rejection, fixed `U/C` rows for indirectly stochastic zero-Q states, time-varying potential-B adjacency rejection when zero-Q `x0` is estimated, and deterministic/indirectly-stochastic state-class stability.

The final package audit also checks maintained Fortran for lines over 132 columns, semicolon-separated statements, legacy real kinds/D-exponents, self-comparison NaN tests, missing dummy `INTENT`/`VALUE` documentation, duplicate Fortran sources, copied dependency source, and build/archive artifacts.
