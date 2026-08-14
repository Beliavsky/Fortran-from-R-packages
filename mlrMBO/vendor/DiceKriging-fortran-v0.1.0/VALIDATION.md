# Validation

## Compiler configuration

Release validation used GNU Fortran 14.2.0 with:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O2
```

All library modules, all permanent test programs, and the example compile and
run successfully with those flags.

The environment used for this translation did not contain the `fpm` executable.
`fpm.toml` was parsed independently as TOML and the exact FPM source/test/example
layout was compiled directly with gfortran.

## Permanent tests

`test_covariance.f90`

- checks DiceKriging Gaussian range scaling against the direct formula;
- checks symmetry and nugget handling;
- checks power-exponential parameter flattening and bounds;
- checks isotropic parameter propagation and scaling-function derivatives.

`test_upstream_regression.f90`

- reproduces the upstream 16-point Branin Matern-5/2 constant-trend regression:
  ranges approximately `(0.8254355, 2.0)`, variance `145556.6`, trend
  coefficient `306.5783`;
- reproduces the upstream interaction-trend regression:
  ranges approximately `(0.7917705, 2.0)`, variance `87350.78`, coefficients
  `(579.5111, -402.8916, -362.0008, 431.2314)`;
- compares all 16 LOO means and standard deviations with upstream test targets.

`test_gradients.f90`

- checks the noisy Gaussian log-likelihood analytic gradient against centered
  finite differences;
- checks the LOO analytic gradient against centered finite differences.

`test_pmle_scaling.f90`

- reproduces the upstream known-nugget SCAD/PMLE fit (`beta=-0.5586176`,
  `sd2=3.35796`, range `2.417813`, nugget `0.001`);
- reproduces the upstream nonlinear-scaling fit with eta values approximately
  `(17.6113829, 2.4169448, 0.8873958)`.

`test_update_prediction.f90`

- checks interpolation at design points;
- checks 95% prediction-interval construction;
- verifies that `cov_reestimate=.true.` changes the fitted range after an
  appended observation;
- verifies response-only update behavior.

## Additional development checks

Before the permanent suite was consolidated, the following were also checked
under runtime bounds/allocation checking:

- stationary Matern-5/2 likelihood gradient versus finite differences;
- estimated-nugget profile-likelihood gradient versus finite differences;
- noisy likelihood gradient and LOO gradient separately;
- Hartman-3 and Hartman-6 values against an independent implementation;
- exact upstream Branin LOO regression arrays;
- clean optimized rebuild from an empty build directory.

The strongest compatibility check is the upstream Branin MLE regression: the
Fortran implementation obtains ranges `0.8254355` and `2.0000000`, process
variance about `145556.6`, and trend coefficient `306.5783`, demonstrating that
the covariance parameterization, profile likelihood, regression profiling,
bounds, and optimizer agree on the same fitted solution.
