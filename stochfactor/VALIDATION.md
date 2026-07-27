# Validation

## Environment

- GNU Fortran 14.2.0
- GNU `make`
- LAPACK and BLAS from the system linker
- No R installation
- No `fpm` installation

## Debug build

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -O0 -g -fcheck=all -fbacktrace
```

Command:

```sh
make debug
```

Validated:

- Compilation of every library module.
- Univariate SV simulation, stationary initialization, likelihood, Omori mixture,
  MCMC, regression, Student-t, leverage, residual, prediction, and rolling paths.
- Factor-SV simulation, zero-factor handling, covariance/correlation paths,
  initialization, restrictions, fitting, Normal-Gamma shrinkage, identification,
  Woodbury calculations, prediction, likelihoods, and densities.
- Demo, CSV, and custom-restriction applications.
- License headers in every Fortran source/test/application file.

Result:

```text
Univariate SV simulation, likelihood, MCMC, regression, t-error, leverage, residual, and prediction tests passed.
Factor-SV simulation, covariance, initialization, fitting, identification, prediction, and density tests passed.
GPL-2.0-or-later source license checks passed.
```

## Optimized build

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -O2 -fbacktrace
```

Command:

```sh
make release
```

The same tests and applications passed.

## Numerical checks of particular importance

- Omori mixture probabilities sum to one.
- Simulated conditional volatilities and Student-t scales are positive.
- MCMC parameter arrays remain finite and constrained (`sigma > 0`, `nu > 2`,
  `abs(rho) < 1`).
- Regression and rolling-forecast paths return finite values.
- Factor covariance matrices are symmetric and factor correlation diagonals equal
  one within tolerance.
- Custom loading restrictions remain exactly zero.
- Zero-factor simulation and fitting execute without indexing failures.
- Normal-Gamma local and global shrinkage variables stay positive.
- Woodbury precision equals the precision stored by the prediction routine.
- Covariance-form and Woodbury-form predictive log likelihoods agree within
  `1e-8` in the regression test.
- Multivariate Normal log density is checked against the identity-covariance
  analytic value.

## Equivalence statement

R was unavailable, so no direct R-versus-Fortran chain comparison was run. MCMC
chains are inherently sensitive to RNG, parameterization, and transition kernels.
The tests establish internal formula consistency, constraints, array safety, and
executable end-to-end workflows; they do not establish exact posterior-draw
identity with the original packages.
