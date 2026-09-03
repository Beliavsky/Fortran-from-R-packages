# Build and validation record

Validation date: 2026-09-02

## Environment

- GNU Fortran: `GNU Fortran (Debian 14.2.0-19) 14.2.0`
- FPM: not installed in this execution environment
- fprettify: not installed in this execution environment
- Windows cross-gfortran: not installed in this execution environment

The environment also could not download missing build tools, so it was not
possible to install FPM or fprettify during this run.

## Compilation actually performed

A fresh temporary build directory was created with no pre-existing `.o` or
`.mod` files.  The complete sibling TMB source set and all maintained glmmTMB
modules were compiled with:

```text
-std=f2018
-Wall
-Wextra
-Wimplicit-interface
-Werror=implicit-interface
-fcheck=all
-pedantic
```

The deterministic test program and example were then linked from those freshly
compiled objects and executed.

Result:

```text
All glmmTMB deterministic tests passed.
NB2 observation log-likelihood:    -2.185170
Random-effect negative log-likelihood:     1.755310
Implied correlation:     0.287348
```

No `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent finite-only
assumption was used.

Compiler diagnostics consisted of intentional real-equality warnings for exact
count/boundary checks and pre-existing warnings in the sibling TMB translation.
There were no implicit-interface errors and no runtime-check failures.

## Deterministic coverage exercised

The test suite checks:

- robust beta-binomial and negative-binomial log masses;
- generalized-Poisson and skew-normal densities;
- Bell/Lambert-W distribution support;
- COM-Poisson Poisson-limit density and variance;
- positive and zero-mass Tweedie likelihoods;
- arbitrary-smoothness Matérn correlation against an independent value;
- one-dimensional Wishart and inverse-Wishart reference densities;
- LKJ transformed-parameter kernel;
- stable link transformations;
- zero-inflated Poisson likelihood for zero and positive outcomes;
- family variance helpers;
- unstructured, compound-symmetry, AR(1), OU, and exponential-spatial random
  effect likelihoods;
- reduced-rank loading construction and spherical-to-data-scale transform;
- normal, gamma, Student-t, and Cauchy prior kernels;
- dense linear-predictor construction;
- composition of observation and random-effect terms into the joint NLL.

## Mechanical source/package audits

The final maintained source tree passed checks for:

- free-form lines no longer than 132 columns;
- no semicolon-separated Fortran statements;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no self-comparison NaN tests;
- explicit `INTENT` or `VALUE` on every dummy argument;
- one dummy argument declaration per declaration line;
- a trailing meaningful `!!` FORD comment on every dummy declaration;
- no duplicate Fortran source files;
- no copied TMB/R/BLAS/LAPACK dependency source;
- no `-llapack`, `-lblas`, or system BLAS/LAPACK FPM link directives;
- no new `rfortran-compat` dependency;
- no object, module, static-library, executable, cache, or ZIP products inside
  the package directory;
- successful parsing of `fpm.toml` as TOML.

## FPM limitation

The requested literal commands

```text
fpm build
fpm test
fpm clean --all
```

could not be executed because FPM is not installed in this environment.  This
is an environment limitation and is not represented as a successful FPM run.
The FPM manifest uses only the sibling path dependency:

```toml
TMB = { path = "../TMB" }
```

The corresponding complete source, test, and example targets were instead
compiled and executed directly with GNU Fortran as described above.

Because no Windows compiler/cross-compiler is installed here, Windows execution
was not directly tested.  Maintained sources use standard free-form Fortran and
no POSIX APIs, and the manifest contains no system-library linker directives.
