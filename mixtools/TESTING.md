# Testing

## FPM

```text
fpm test
fpm run
```

## Direct GNU Fortran validation

On Unix-like systems:

```text
./run_gfortran_tests.sh
```

On Windows:

```text
run_gfortran_tests.bat
```

The strict configuration uses Fortran 2018 conformance, warnings as errors,
runtime bounds checking, floating-point traps, and backtraces. Release
validation uses `-O3 -Werror`.

## Test programs

- `test_distributions`: densities, weighted statistics, matrix square roots,
  ellipses, and permutations
- `test_parametric`: normal, gamma, multivariate-normal, multinomial, and
  repeated-measures mixtures
- `test_regression`: linear, logistic, Poisson, expert-gated, and grouped
  regression mixtures
- `test_semiparametric_reliability`: product-kernel EM, symmetric mixtures,
  censored reliability models, model selection, bootstrap, and MCMC

Tests use deterministic synthetic data and the package's portable RNG.

## Cross-validation against the R package

### Multivariate-normal mixture, 2026-08-23

The comparison workflow generated 600 two-dimensional observations from a
two-component normal mixture with true weights `0.4` and `0.6`. The R package
and Fortran translation fitted the same observations using equivalent initial
weights, means, and pooled covariance matrices. Components were ordered by
their first fitted mean before comparison.

Run the complete workflow with either:

```text
run_mvnormal_comparison.bat
bash run_mvnormal_comparison.sh
```

The validation run used:

- R 4.6.1 and CRAN `mixtools` 2.0.0
- FPM 0.12.0 alpha
- GNU Fortran 17.0.0 experimental, 2026-05-10 build
- Windows, with the Fortran executable already compiled

Numerical results were:

| Quantity | R | Fortran | Absolute difference |
| --- | ---: | ---: | ---: |
| Log-likelihood | -1982.23050557 | -1982.230505951 | 3.81e-7 |

Across the fitted parameters, the maximum absolute differences were `3.60e-6`
for component weights, `2.12e-5` for means, and `4.33e-5` for covariance-matrix
entries.

Elapsed wall-clock timings from that run were:

| Stage | R (seconds) | Fortran (seconds) | R/Fortran ratio |
| --- | ---: | ---: | ---: |
| Read observations | 0.000649 | 0.004492 | 0.14 |
| Fit mixture | 2.695228 | 0.024856 | 108.4 |
| Overall startup-through-fit | 4.434529 | 0.029349 | 151.1 |

The fit time is the most relevant performance comparison. R's overall time
includes loading the `mixtools` package, while neither overall measurement
includes process-launch time and the Fortran measurement excludes compilation.
The timings are from one run on one small dataset and are not a controlled or
general-purpose benchmark. Repeated warm runs, larger datasets, additional
dimensions and component counts, difficult initializations, and nearly
singular covariance matrices should be tested before drawing broad performance
conclusions.

This comparison exposed and led to correction of an `mvnormalmixEM` convergence
defect: the translated E-step overwrote the previous log-likelihood, causing
the EM loop to stop after one iteration. `test_parametric` now requires the
multivariate-normal EM fit to perform more than one iteration, preventing that
specific regression.
