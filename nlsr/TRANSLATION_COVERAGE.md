# Translation coverage

## Translated computational code

- `nlfb` stabilized nonlinear least-squares iteration.
- Adaptive Marquardt damping (`lamda`, `laminc`, `lamdec`).
- `phi` and `psi` augmented-Jacobian stabilization.
- QR least-squares step and upstream-style residual orthogonality (`roff`) test.
- Backtracking with `stepredn` / `nbtlim`.
- Lower/upper bounds and exact masks (`lower == upper`).
- Active-bound gradient logic (`bdmsk` values `-3`, `-1`, `0`, `1`).
- Fixed and callback-generated residual weights.
- Forward, backward, central, and Richardson residual Jacobians.
- `resgr`/`resss` numerical roles.
- Covariance / standard-error calculation from the weighted Jacobian.
- The numeric model, Jacobian, and initializer from `SSlogisJN`.

## Represented differently

### `nlxb`, `model2rjfun`, `model2ssgrfun`

These routines mainly translate R formulas/environments into residual/Jacobian functions.  Fortran has no R formula language, so the Fortran API directly accepts typed residual and Jacobian procedures.  This preserves their numerical role without embedding an R parser.

### Symbolic derivative system

`dex`, `nlsDeriv`, `codeDeriv`, `fnDeriv`, `newDeriv`, `newSimplification`, `nlsSimplify`, and `findSubexprs` manipulate R parse trees and environments.  They are R-language metaprogramming rather than a numerical algorithm.  Analytic Fortran callbacks provide the direct equivalent; numerical Jacobian routines cover models without an analytic derivative.

### `wrapnlsr` / `nlsr`

Upstream `wrapnlsr` runs `nlxb` and then calls R's `stats::nls`.  The external second-stage R solver is not part of nlsr's numerical implementation and is not reproduced.  The translated `nlfb` solver is standalone.

### Dynamic weights

The Fortran callback is evaluated using the current parameter vector and current raw residual.  Upstream's closure retains a `resraw` variable whose update behavior is tied to R scoping; the Fortran definition is explicit and deterministic.

## Omitted as non-computational R infrastructure

- S3 `print`, `predict`, `fitted`, `coef`, and model-object methods.
- Formula parsing and environment lookup.
- `digest`-based expression bookkeeping.
- Console tracing/formatting helpers.
- `nlsrSS`'s R `selfStart`/`getInitial` orchestration; the actual logistic initializer is translated.
