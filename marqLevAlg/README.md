# marqLevAlg-fortran

Modern Fortran/FPM translation of the computational core of the R package
**marqLevAlg 2.0.8**.

## Implemented

- Marquardt-Levenberg/Newton optimization with adaptive diagonal inflation.
- Minimization and maximization.
- Numerical objective derivatives (`deriva` semantics).
- Analytical gradient with numerical information matrix (`deriva_grad`).
- User-supplied information/Hessian matrix path.
- Logarithmic one-dimensional step search (`searpas`).
- Parameter, objective, and relative-distance-to-minimum/maximum (RDM)
  convergence tests.
- Optional partial-Hessian RDM convergence.
- `multipleTry` starting-point halving for non-finite initial numerical
  derivatives.
- `blinding` handling of non-finite trial objective values.
- Linear mixed-model random-intercept log-likelihood and analytical gradient
  (`loglikLMM`, `gradLMM`).
- Modern packed-symmetric matrix helpers and SPD solve/inversion routines.

## Build

```sh
fpm build
fpm test
```

## Basic use

```fortran
use marqlevalg, only : dp, mla_result, marqlev_optimize

type(mla_result) :: fit
real(dp) :: x0(2)

x0 = [8.0_dp, 9.0_dp]
call marqlev_optimize(x0, objective, fit)
```

Overloads are available for objective-only, objective+gradient, and
objective+gradient+information/Hessian callbacks.

## Important Hessian convention

The upstream R source copies a user-supplied `hess()` result directly into the
matrix slot that its Marquardt step treats as the **information matrix**.  This
translation preserves that behavior.  For a mathematically conventional
Hessian `H` of the user's objective, the matrix required by the internal
maximization algorithm is `-H_internal`; for ordinary minimization this is
`+H_user`.

## Parallelism

The R package uses `foreach`/`doParallel` only to evaluate independent finite-
difference points concurrently.  The Fortran library is serial and standalone;
the derivative formulas and optimizer are unchanged.  No R runtime is needed.

See `TRANSLATION_COVERAGE.md` for exact differences and `VALIDATION.md` for the
test configuration.
