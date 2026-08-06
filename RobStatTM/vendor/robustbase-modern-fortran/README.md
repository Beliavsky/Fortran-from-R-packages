# robustbase modern Fortran

A modern Fortran numerical translation and reimplementation of the computational core of the R package `robustbase` 0.99-7.

The project uses plain arrays, procedures, and derived result types. Plotting, R formulas, model frames, S3/S4 classes, factors/contrasts, and R time-series or missing-value metadata are intentionally excluded.

This release substantially expands the earlier translation. It does not claim byte-for-byte equivalence with the R, C, or legacy Fortran implementation. Every feature listed under **Implemented and tested** is present in the source and exercised by the test suite.

## Requirements

- A Fortran 2018 compiler
- LAPACK
- BLAS

Validation used GNU Fortran 14.2.0.

## Build and test

```sh
make check
make release-check
```

Equivalent direct commands are:

```sh
./test/run_all.sh debug
./test/run_all.sh release
```

The debug build uses runtime bounds and consistency checks. Both configurations treat compiler warnings as errors.

An `fpm.toml` manifest is included and links LAPACK and BLAS. `fpm` was unavailable in the validation environment, so that build route is not claimed as tested.

## Implemented and tested

### Robust scales and location

- Exact-definition `Qn` and `Sn`, including default finite-sample corrections
- Weighted high median
- MAD and IQR scales
- Huber location and scale outputs
- Huberization
- Huber tau and tau-style scales
- Trimmed mean absolute deviation

### Score, loss, and weight functions

- Huber
- Hampel
- Tukey bisquare
- Welsh/Gaussian weight
- Optimal redescending score
- Generalized Gaussian weight (GGW)
- Linear-quadratic-quadratic (LQQ)
- First derivatives for all implemented score families
- LQQ second derivative
- Numerical normalized GGW loss integration
- Hard-rejection and smooth weights

### Medcouple and adjusted boxplots

- Medcouple with tie handling
- Left and right medcouple
- Adjusted boxplot fences and outlier flags

### Matrix utilities

- Row and column medians
- Robust centering and scaling
- Independent-column detection
- Full-rank matrix reduction
- SVD rank using the `rankMM` tolerance rule
- Classical PCA with optional centering, scaling, scores, rank, and deterministic loading signs

### Robust covariance and outlyingness

- Comedian and comedian covariance
- Iterative comedian covariance
- Gnanadesikan-Kettenring covariance
- Orthogonalized Gnanadesikan-Kettenring covariance
- Random-start concentration-step MCD
- Deterministic six-start MCD
- Singular exact-fit hyperplane output from deterministic MCD
- Partitioned FAST-style MCD with full-data concentration refinement
- Analytical raw and reweighted MCD consistency and finite-sample factors
- Robust Mahalanobis distances
- Projection outlyingness
- Full continuous-data adjusted outlyingness with sampled hyperplane directions
- Tolerance-ellipse coordinates

### Robust linear regression

- Basic random-start LTS
- Advanced LTS with exhaustive, deterministic, random, and automatic search modes
- Partitioned FAST-style LTS with full-data refinement
- Raw and reweighted LTS estimates, scales, weights, covariance, objective, and subset
- LTS-initialized MM regression
- Least-absolute-residual (`lmrob.lar`) numerical fit
- S regression with simple, nonsingular, exhaustive, and automatic subsampling
- SM/MM chains
- SMDM chains
- Weighted, AVAR1-style, and sandwich covariance estimates
- Robust residual weights, scales, fitted values, standard errors, and exact-fit flags

### Robust generalized linear models

- Huberized Pearson IRLS for binomial and Poisson models
- Bianco-Yohai robust logistic regression and analytical derivatives
- Mqle-style binomial and Poisson estimation with bias correction and sandwich covariance
- MT-style binomial and grouped-binomial estimation using variance-stabilizing transformations

### Robust nonlinear regression

- Generic finite-difference robust nonlinear least squares
- MM estimation
- Tau estimation
- Constrained-M estimation
- Maximum-trimmed-likelihood estimation
- Bounded optimization, robust covariance matrices, fitted values, residuals, and weights

### Inference and diagnostics

- Robust confidence and prediction intervals
- Wald tests for arbitrary coefficient subsets
- Robust deviance comparisons with user-specified degrees of freedom
- Robust R-squared
- Leverage and standardized-residual outlier statistics
- Numerical Hessian/covariance outputs where applicable

## CSV application

The application reads a headerless or headed two-column file containing predictor and response values.

```sh
build/debug/bin/fit_csv data/example_xy.csv lts
build/debug/bin/fit_csv data/example_xy.csv fastlts
build/debug/bin/fit_csv data/example_xy.csv partlts
build/debug/bin/fit_csv data/example_xy.csv lar
build/debug/bin/fit_csv data/example_xy.csv lmrob
build/debug/bin/fit_csv data/example_xy.csv smdm
build/debug/bin/fit_csv data/example_binary.csv by
build/debug/bin/fit_csv data/example_binary.csv mqle-binomial
build/debug/bin/fit_csv data/example_binary.csv mt
```

Poisson Mqle mode is also available as `mqle-poisson` for nonnegative count responses.

## Numerical differences from robustbase

- `Qn`, `Sn`, and the medcouple use direct exact-definition arrays rather than the package's optimized selection/search kernels. Their asymptotic memory and runtime are therefore less favorable.
- Partitioned MCD and LTS preserve concentration-step and partition/refinement ideas but are clean numerical reimplementations, not line-for-line ports of every optimized kernel and control heuristic.
- Historical simulation-table overrides for a few MCD finite-sample correction combinations are not embedded. The analytical correction formulas are implemented.
- `lmrob` S/M/D paths are plain-array reimplementations. Exact-design splitting for categorical model matrices, every leverage correction, and the complete expert-control surface are not reproduced.
- LAD/LAR uses a smoothed iteratively reweighted solver and a quantile-regression covariance approximation rather than an external linear-programming solver.
- `glmrob` and `nlrob` methods preserve their principal estimating criteria and workflows but use self-contained numerical solvers. Exact R optimizer endpoints are not expected.
- GGW rho is evaluated by deterministic Simpson integration for general tuning constants rather than R's adaptive integration routine or the original hard-coded polynomial tables.
- Intrinsic Fortran random streams replace R random streams.

## Remaining exclusions

The remaining exclusions are primarily non-numerical or exact-compatibility surfaces:

- R formulas, model frames, offsets, factors, contrasts, family objects, and missing-value policies
- S3/S4 objects and methods, summaries, printing, plotting, and package-data wrappers
- Categorical-data handling in adjusted outlyingness
- Exact historical correction lookup tables and every legacy warning/control code
- Automatic calibration routines that solve for GGW/LQQ tuning constants from requested efficiency or breakdown point
- Exact iteration-by-iteration reproduction of legacy C/Fortran optimized FAST-MCD, FAST-LTS, `lmrob`, `glmrob`, and `nlrob` kernels

See `API_MAP.md` for a routine-level mapping and `VALIDATION.md` for the executed test scope.

## License

The original package declares `GPL (>= 2)`. This translation is distributed under **GPL-2.0-or-later**. `LICENSE` contains GNU GPL version 2, and every Fortran source, application, example, and test file carries the corresponding SPDX identifier and license notice.
