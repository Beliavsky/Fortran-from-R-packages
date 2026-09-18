# Validation record

Validation performed while preparing this translation and compatibility update:

- Compiler: GNU Fortran 14.2.0.
- Maintained Fortran sources were compiled from a clean temporary directory with:

  ```text
  -std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wconversion-extra -Werror=line-truncation -pedantic
  ```

- The package sources compiled successfully against local API-compatible validation shims for the exact shared procedures imported from `rfortran-core`, `rfortran-linalg`, and `MTS`. Those shims are validation-only and are not present in this package.
- `test/test_steadystatebvar.f90` compiled and ran successfully end-to-end.
- `example/basic_workflow.f90` compiled and ran successfully.
- `example/stochastic_volatility.f90` compiled and ran successfully.

## Deterministic test coverage

The test suite exercises:

- `bvar_create` and all three deterministic setup variants;
- OLS setup quantities and future deterministic-term construction;
- Minnesota and steady-state priors, positivity validation, and beta restrictions;
- homoscedastic Jeffreys and inverse-Wishart posterior branches;
- fixed-seed repeatability for the native homoscedastic sampler and predictive recursion;
- unconditional forecasts and annualized-growth transformations;
- equality-conditioned homoscedastic forecasts;
- homoscedastic OIRFs and GIRFs;
- RW stochastic-volatility prior validation, posterior sampling, positive `phi` draws, complete in-sample covariance paths, future covariance draws, forecasts, time-specific IRFs, and summaries;
- AR(1) stochastic-volatility posterior sampling, `gamma_1` stationarity bounds, `Phi` draws, complete covariance paths, forecasts, time-specific GIRFs, and summaries;
- the upstream-compatible rejection of stochastic-volatility conditional forecasts;
- `ppi` and mean/median posterior summaries.

The compatibility update additionally checked the SV equations against the retained upstream RW/AR1 Stan programs: `A` construction, centered latent-state transition densities corresponding to the Stan noncentered parameterization, RW inverse-gamma updates, AR1 inverse-Wishart updates, `Sigma_u,t = A^{-1} Lambda_t A^{-T}`, and the generated predictive volatility recursion.

## Static checks

The final package tree passed checks for:

- free-form maintained Fortran only;
- no maintained source line longer than 132 characters;
- no semicolon-separated Fortran statements;
- no self-comparison NaN tests;
- no `double precision`, `real*8`, or D-exponent literals in maintained Fortran;
- no system BLAS/LAPACK link flags;
- no copied shared-dependency source directories;
- no duplicate maintained Fortran source files;
- no object/module/executable/library/cache/ZIP build products inside the package;
- explicit `INTENT` or `VALUE` on every dummy argument;
- one dummy argument per declaration;
- a meaningful trailing FORD `!!` comment on every dummy-argument declaration;
- 10 translation mapping entries agreeing with `r_functions_translated = 10`, `r_functions_total = 10`, and `frac_functions_translated = 1.0`;
- every `r_source` and `fortran_source` path in `fpm.toml` resolving to a retained file.

## FPM environment limitation

The execution environment used to prepare this artifact does not contain the `fpm` executable. Explicit final attempts to run all three requested commands returned shell status 127 (`fpm: command not found`):

```text
fpm build
fpm test
fpm clean --all
```

Therefore the successful direct GNU Fortran validation above is **not** represented as an FPM run. The package is structured for repository-root placement with sibling dependencies and should receive the final repository-level validation on a machine with FPM installed:

```text
fpm build
fpm test
fpm run --example basic_workflow
fpm run --example stochastic_volatility
fpm clean --all
```
