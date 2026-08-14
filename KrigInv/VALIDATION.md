# Validation

KrigInv-fortran v0.2.0 was rebuilt from a clean directory with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wno-maybe-uninitialized \
-Wimplicit-interface -fcheck=all -O2
```

`-Wno-maybe-uninitialized` suppresses gfortran false positives involving allocatable derived-type results; all other enabled warnings are errors.

## Permanent tests

Nine FPM tests pass:

1. `test_core`
   - interpolation and UK covariance symmetry;
   - excursion probability and bivariate-normal identities;
   - Sobol integration;
   - Ranjan/Bichon/TMSE/TSEE criteria.
2. `test_batch`
   - SUR/IMSE variance reduction;
   - discrete single-point optimization;
   - greedy SUR/TIMSE/Vorob/future-volume batch optimization.
3. `test_reference`
   - independent noisy two-dimensional Matern-5/2 UK mean, standard deviation, and posterior covariance reference values.
4. `test_egi`
   - legacy fixed-parameter sequential EGI and one-step SUR EGI.
5. `test_cov_reestimate`
   - `fit_krig_model` through the vendored DiceKriging backend;
   - published Branin MLE regression target;
   - default covariance re-estimation after appending a point;
   - explicit `cov_reestimate=.false.` preservation of range and variance.
6. `test_egi_reestimate`
   - fitted-model EGI defaults to `CovReEstimate=model@param.estim`;
   - explicit `cov_reestimate=.false.` keeps covariance parameters fixed.
7. `test_noise_reestimate`
   - transition from a deterministic fitted model to a new observation with positive observation variance;
   - noise-vector propagation and covariance refitting;
   - finite UK prediction after the noisy refit.
8. `test_dice_covariance`
   - Gaussian covariance scaling;
   - nugget behavior;
   - power-exponential parameter flattening/bounds;
   - isotropic parameters;
   - nonlinear scaling values and derivatives.
9. `test_dice_gradients`
   - analytic noisy-likelihood gradients versus central finite differences;
   - analytic LOO gradients versus central finite differences.

Both examples also compile and run under the same flags.

## Published DiceKriging regression

The Branin regression used by the translated DiceKriging test suite is exercised through `fit_krig_model`, not by calling `dk_model` directly.  With the same deterministic RNG seed and optimizer settings, the fitted KrigInv model reproduces the upstream DiceKriging regression targets within the test tolerances:

```text
range ~= [0.8254355, 2.0000000]
variance ~= 145556.6
constant trend ~= 306.5783
```

This simultaneously checks the KrigInv-to-DiceKriging adapter, covariance conventions, concentrated likelihood, bounds, GLS trend profiling, and multistart bounded optimizer.

## CovReEstimate regression

For a DiceKriging-backed fitted model:

- `model%param_estim` is true;
- `model%cov_reestimate_default` is initialized from it;
- `update_krig_model` with no explicit covariance flag refits covariance/variance;
- `update_krig_model(...,cov_reestimate=.false.)` keeps them fixed;
- `egi` forwards the same default/override behavior across sequential infill updates.

The noisy update regression additionally confirms that an all-zero default `new_noise` does not unnecessarily switch a noise-free fit into the heteroskedastic-noise likelihood case, while a positive new noise variance does.

## Earlier independent numerical checks retained

The v0.1.0 numerical checks remain applicable to the legacy/model-independent KrigInv calculations:

- the fixed noisy two-dimensional UK case agreed with an independent NumPy implementation to about `4e-16` in mean, standard deviation, and posterior covariance;
- the native bivariate-normal CDF agreed with SciPy over 300 randomized `(a,b,rho)` cases to about `1.4e-16` maximum absolute error.

## FPM

The `fpm` executable is not installed in the validation container.  `fpm.toml` parses successfully as TOML and the exact source/test/example tree was compiled and linked directly with gfortran using the flags above.  FPM's module dependency resolver should compile the same module graph automatically.
