# Validation

## Environment

Validation for this translation was performed with:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Linux x86-64
```

The sandbox did not contain an `fpm` executable, and its network path cannot
resolve GitHub hosts. Consequently the literal commands `fpm build`, `fpm test`,
and `fpm clean --all` could not be executed in this environment. The manifest
was parsed as TOML and is laid out for sibling packages at the repository root.

Current repository sources were checked through the available web source view.
The interfaces used here match `rfortran-core`/`rfortran-linalg`, and the
existing top-level `optimx` package exports `optimx_mod`, `optimx_problem`,
`optimx_control`, `optimx_result`, `initialize_problem`, and `optimr`. Its
callback signature and result/control fields used by `dlm_mle` were checked
against the current repository source. The `optimx` and `rfortran-core` real
kinds both resolve to the gfortran binary64/double kind used by this build.

`rfortran-linalg` uses the repository-standard pinned pure-Fortran LAPACK FPM
dependency, so this package contains no system BLAS/LAPACK link directive.

## Direct GNU Fortran validation

Because the sibling packages were not present locally, the complete maintained
`dlm` source graph was compiled against small API-compatible validation stand-ins
for `r_kinds`, the two `r_linalg` procedures used here, and the public subset of
`optimx_mod` used by `dlm_mle`. These stand-ins live outside the package tree and
are not included in the release archive. The `optimx` stand-in uses a simple
deterministic optimizer only to exercise the DLM objective/builder integration;
the production manifest points to `../optimx`.

The checked build used:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface
-Wconversion-extra -Wcompare-reals -Werror -fcheck=all -fbacktrace
-ffree-line-length-none
```

The optimized build used:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
-Wconversion-extra -Wcompare-reals -Werror -fcheck=all -fbacktrace
-ffree-line-length-none
```

Both configurations compiled all maintained source and produced:

```text
test_arms: PASS
test_bounds: PASS
test_kalman: PASS
test_mle: PASS
test_models: PASS
test_random_mcmc: PASS
```

All examples also compiled and ran in both configurations:

```text
arms_sampling: sample mean = -0.00683
local_level: filtered final state = 1.34350
local_level: smoothed initial state = 0.97726
mle_variance: fitted log variance = 0.69315
```

The ARMS tests cover deterministic seeding, bounded univariate standard-normal
sampling, a deliberately non-log-concave symmetric two-mode target to exercise
the Metropolis correction, and the multivariate random-direction hit-and-run
wrapper. The MLE test uses an IID zero-mean Gaussian DLM whose analytic optimum
is `log(2)` and checks the recovered parameter and objective value. Existing
tests continue to cover constructors, stationary AR transformation, time-varying
matrix indexing, missing-observation filtering, likelihood, smoothing,
forecasting, convex bounds, MCMC summaries, Wishart simulation, seed
reproducibility, and `dlmGibbsDIG`.

No GNU executable-stack warning remains in the final MLE implementation; its
optimizer callback uses a synchronous module context rather than an internal
procedure trampoline. Reentrant/concurrent `dlm_mle` calls are rejected by this
bridge.

## Source audit

Before packaging, all maintained Fortran under `src/`, `test/`, and `example/`
was checked. The audit found:

- only free-form `.f90` maintained Fortran source;
- ASCII maintained source and no maintained line longer than 132 characters;
- a single package real kind `dp` imported from `r_kinds`;
- no `double precision`, `real*8`, `kind(0.0d0)`, or `d0`/`D0` exponent forms;
- no semicolon-separated executable statements;
- no self-comparison NaN tests; missing-data checks use `ieee_is_nan`;
- no fast-math or finite-math-only flags;
- no system BLAS/LAPACK link directives;
- no duplicate maintained Fortran source files;
- every dummy argument declared separately with `INTENT` or `VALUE` and a
  trailing `!!` FORD documentation comment;
- 25 mapping entries for 25 exported computational R functions, with an empty
  untranslated list in `fpm.toml`;
- no copied `optimx`, `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, ARPACK,
  or other dependency source inside the package;
- no object files, module files, executables, caches, build directories, or ZIP
  files inside the package tree.

`fprettify` was not installed in the environment, so it could not be run. The
new source was formatted manually to the same free-form style and audited for
line length and the requested declaration rules.

Because FPM was unavailable, the requested `fpm clean --all` step was replaced
by an explicit package-tree artifact scan immediately before archive creation.
No FPM build directory was created inside the package.
