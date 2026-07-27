# Validation status

## Completed smoke checks

The project was compiled with GNU Fortran 14.2 using:

```console
gfortran -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
```

The checked test suite currently covers:

- Student-t, skew-GED, and Johnson-SU CDF/quantile round trips.
- Scaled modified-Bessel K values against the exact order-one-half formula.
- NIG equivalence to generalized hyperbolic `lambda = -1/2`.
- Symmetric GH skew-t equivalence to standardized Student-t.
- GHYP, NIG, and GH skew-t CDF/quantile round trips and random generation.
- GJR-GARCH simulation, filtering, likelihood, fitting, and forecasts.
- Generalized FIGARCH polynomial weights and FIGARCH(2,d,1) fitting.
- Component-GARCH simulation and fitting.
- realGARCH simulation, filtering, joint likelihood, fitting, and measurement
  residual output.
- Hentschel ALLGARCH simulation/filtering/fitting and construction of all eight
  supported fGARCH submodels.
- Standard GARCH fitting on generated CSV data.
- Fractional differencing/integration checks and ARFIMA simulation/fitting.
- VaR coverage and directional-accuracy test execution.
- The no-argument demonstration and CSV fitting example.

## Interpretation

These are deterministic identities, generated-data smoke checks, and runtime
safety tests. They establish that the implemented paths compile and execute and
that several internal mathematical identities hold. They are not a complete
cross-language validation suite.

## Not established

- Full numerical equivalence to R `rugarch` 1.5-6 across all models and data.
- Equivalence to R's solver choices, initialization, scaling, constraints,
  gradients, Hessians, or standard errors.
- Accuracy of every GH-family tail probability at extreme parameter values.
- Correctness for every boundary case or non-Gaussian model combination.
- Cross-compiler or cross-platform behavior beyond GNU Fortran in this build.
- Production suitability.

The test suite should be treated as a smoke suite, not a certification suite.


## Version 0.3 smoke validation

The strict suite additionally exercises:

- weighted portmanteau and ARCH-LM diagnostics;
- Nyblom, sign-bias, and adjusted Pearson statistics;
- Newey-West and inverse-Hessian covariance calculations;
- stationary and fixed-block bootstrap indices;
- Weibull VaR-duration, GMM, Hong-Li, and MCS tests;
- multi-series fitting and forecasting;
- rolling one-step forecasts;
- parametric simulation/refit distributions; and
- ARFIMA forecasting and automatic order selection.

These are implementation smoke tests, not exhaustive numerical equivalence
checks against every R `rugarch` option.


## Version 0.4 smoke validation

The additional checked test exercises:

- external mean and variance regressors;
- variance targeting and extended filtering;
- numerical Hessian and Newey-West covariance attachment;
- partial raw/kernel/semi-parametric and full simulation/refit bootstrap paths;
- serial ARFIMA rolling, bootstrap, parameter-distribution, multi-series, and
  cross-validation workflows;
- forecast performance, GH parameter transformation, distribution moments, and
  numeric forward-index utilities.

These remain smoke checks. The two-step external-regressor estimator, numerical
covariance, and dependency-free semi-parametric bootstrap are not claimed to be
bit-for-bit equivalents of the R solver and its optional package dependencies.
