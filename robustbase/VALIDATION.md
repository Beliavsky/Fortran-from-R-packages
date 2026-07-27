# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018
- LAPACK and BLAS
- Linux x86-64

## Strict configurations

Debug:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace -ffree-line-length-none
```

Release:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
-Werror -fbacktrace -ffree-line-length-none
```

## Commands

```sh
./test/run_all.sh debug
./test/run_all.sh release
```

The license audit verifies `GPL-2.0-or-later` SPDX and notice text in every `.f90` file under `src`, `app`, `example`, and `test`, and verifies the project license metadata.

## Numerical coverage

### Scale and score suite

- `Qn`, `Sn`, MAD, weighted high median, Huber location, medcouple, adjusted boxplots
- Huber, Hampel, Tukey, Welsh, optimal, GGW, and LQQ score functions
- Finite-difference derivative checks for Welsh, optimal, GGW, and LQQ
- GGW and LQQ loss saturation and score/weight identities
- PCA/rank and matrix utility paths

### Covariance suite

- Comedian, iterative comedian, GK, OGK, random-start MCD
- Deterministic six-start MCD
- Singular exact-fit hyperplane result
- Partitioned FAST-style MCD and correction formulas
- Robust Mahalanobis and adjusted outlyingness

### Regression suite

- Basic and advanced LTS
- Partitioned FAST-style LTS
- MM and least-absolute-residual regression
- S and SMDM `lmrob` chains
- Robust prediction, Wald/deviance comparisons, R-squared, and outlier statistics
- Bianco-Yohai logistic regression
- Mqle binomial and Poisson fits
- MT binomial fit
- Robust nonlinear least squares plus MM, tau, CM, and MTL fits

### Applications

The test driver executes:

- `demo_robustbase`
- `nonlinear_example`
- `next_batch_example`
- CSV modes `mm`, `lts`, `fastlts`, `partlts`, `lar`, `lmrob`, `smdm`, `by`, `mqle-binomial`, and `mt`

## Acceptance criteria

- No compiler warnings in either configuration
- No runtime-check failures in the debug configuration
- Finite covariance/standard-error outputs where claimed
- Robust estimates remain close to known simulated coefficients under contamination
- Exact count/subset/reconstruction identities where appropriate
- Symmetric covariance matrices and closed tolerance ellipses
- All applications exit successfully

## Equivalence limits

R was unavailable in the validation environment. Tests therefore use mathematical identities, hand-derived values, deterministic regression cases, contamination-resistance checks, and cross-method invariants. Exact R random streams, optimizer endpoints, and iteration-by-iteration legacy-kernel equivalence are not claimed.
