# Porting notes

## Scope

The R package is primarily a wrapper around the native NLopt library. This port
translates the R-level computational API to typed Fortran and supplies a
portable self-contained solver backend. R list/S3 printing, dynamic-language
argument forwarding, and the native C registration layer are not reproduced.

## Callback convention

Objective callbacks have the form

```fortran
subroutine objective(x, value, gradient, need_gradient, status)
```

Constraint callbacks similarly return a vector and Jacobian. If an analytic
derivative is unavailable, leave the `intent(inout)` derivative array unchanged
when `need_gradient` or `need_jacobian` is true. The library initializes it to
IEEE NaNs and automatically falls back to central finite differences.

## Algorithm mapping

| Public name | Fortran backend |
|---|---|
| `lbfgs` | projected full-memory BFGS |
| `varmetric` | projected BFGS |
| `tnewton` | projected BFGS |
| `slsqp` | BFGS with outer quadratic penalties |
| `mma` | BFGS with outer quadratic penalties |
| `ccsaq` | BFGS with outer quadratic penalties |
| gradient `auglag` | BFGS with outer quadratic penalties |
| `neldermead` | Nelder-Mead |
| `sbplx` | common Nelder-Mead implementation |
| `bobyqa` | common Nelder-Mead implementation |
| `newuoa` | common Nelder-Mead implementation |
| `cobyla` | coordinate pattern search with penalties |
| derivative-free `auglag` | pattern search with penalties |
| `direct`, `direct_l` | Halton sampling plus local refinement |
| `crs2lm`, `isres` | Halton sampling plus local refinement |
| `stogo`, `mlsl` | Halton sampling plus local BFGS refinement |

These methods preserve the optimization model and high-level calling pattern,
but do not claim to reproduce NLopt's specialized C/C++ algorithms, iteration
counts, random populations, or exact stopping trajectories.

## Bounds

Use IEEE infinities when convenient. Very large finite values such as
`huge(1.0_dp)` are also accepted as practical unbounded limits by local methods.
Global methods require genuinely finite, reasonably scaled bounds.

## Status values

The standard NLopt integer status convention is retained:

- positive values indicate successful termination or a stopping condition;
- `-2` indicates invalid arguments;
- `-4` indicates roundoff-limited progress;
- `-5` indicates a callback-requested stop or nonfinite callback value.

## Licensing

The Fortran translation is LGPL-3.0-or-later, matching the R package. The
original R wrapper sources are retained under `original/`. The bundled NLopt C
source archive is not incorporated into the compiled Fortran library.
