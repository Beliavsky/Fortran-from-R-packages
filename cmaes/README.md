# cmaes-fortran

Modern Fortran translation of the computational code in the R package
**cmaes 1.0-12**, a covariance matrix adapting evolution strategy (CMA-ES)
implementation by Heike Trautmann, Olaf Mersmann, and David Arnu.

The original package states `License: GPL-2`. This translation is distributed
under GPL-2.0-only. See `COPYING`, `LICENSES.md`, and `original/`.

## Implemented computational functionality

- Hansen-style CMA-ES generational update used by `cma_es()`
- covariance matrix adaptation, evolution paths, and global step-size adaptation
- default `lambda`, `mu`, recombination weights, `mueff`, `ccum`, `cs`,
  covariance learning rate, and damping formulas from the R source
- finite and infinite box bounds
- the package's squared-distance multiplicative constraint penalty
- `fnscale`, including maximization for negative `fnscale`
- `stopfitness`, `maxit`, the undocumented `stop.tolx`, and flat-fitness escape
- `keep.best`
- scalar objective callbacks and package-style vectorized population evaluation
- diagnostic histories for sigma, covariance eigenvalues, selected populations,
  and population fitness values
- `extract_population`
- benchmark/helper computations: `f_sphere`, `f_rand`, `f_rosenbrock`,
  `f_rastrigin`, shifted, rotated, and biased objective evaluation

LAPACK `DSYEV` replaces R's `eigen(C, symmetric=TRUE)`. A small deterministic
Park-Miller/Box-Muller RNG replaces R's global RNG. Therefore a numeric seed is
reproducible within this Fortran library but does not reproduce R's random
stream for the same seed.

## Basic use

```fortran
use cmaes, only : dp, cma_control, cma_result, cma_es
use cmaes_functions, only : f_sphere

type(cma_control) :: control
type(cma_result) :: result
real(dp) :: par(3), lower(3), upper(3)

par = [3.0_dp, -2.0_dp, 4.0_dp]
lower = -10.0_dp
upper = 10.0_dp
control%seed = 12345
control%maxit = 500
control%stopfitness = 1.0e-10_dp

result = cma_es(par, f_sphere, lower, upper, control)
```

For a vectorized objective, set `control%vectorized=.true.` and pass the
optional population callback. Its input has shape `(dimension, lambda)` and it
returns `lambda` values.

## Package fidelity notes

The translation follows the executable R source, including several details
that are easy to change accidentally:

- `sigma` is scalar. Although the R documentation mentions a vector, the
  source explicitly checks `length(sigma) == 1`.
- Trial points are clamped before objective evaluation, but covariance and mean
  updates use the original, unclamped trial coordinates, exactly as in R.
- Constraint penalty is `1 + sum((arx-vx)^2)`. A point is considered fully
  valid when that floating-point quantity compares `<= 1`.
- The `stop.tolx` test uses `sigma*pc < tolx` without `abs(pc)`, matching the
  source.
- Flat-fitness detection uses exact floating-point equality.
- Diagnostic `value` history preserves the R source's second indexing by
  `aripop` after sorting `arfitness`; it is not silently rewritten as the
  first `mu` sorted values.
- The package's `f_rosenbrock` first forms `z=x+1`, so its optimum is at
  `x=0`, not at the conventional all-ones point.
- `cmaES()` is only a deprecated R naming wrapper. It is omitted because a
  Fortran module named `cmaes` cannot also contain a procedure of the same
  case-insensitive name.

R parameter names, R list/S3 result objects, warnings/messages, documentation
printing, and other language-interface machinery are not computational and are
not reproduced.

## Build

With FPM and system BLAS/LAPACK:

```text
fpm build
fpm test
fpm run --example sphere
```

The manifest links `lapack` and `blas`.

## Validation

The regression tests cover:

1. stop-fitness convergence on a sphere problem;
2. maximization through `fnscale=-1` and a boundary optimum;
3. all diagnostic histories and population extraction;
4. scalar versus vectorized trajectory equivalence;
5. benchmark/helper functions and deterministic `f_rand` seeding.

See `TRANSLATION_NOTES.md` for a source-level mapping.
