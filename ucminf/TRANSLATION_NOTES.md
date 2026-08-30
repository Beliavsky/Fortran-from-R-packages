# Translation notes

## Included computational behavior

- UCMINF inverse-Hessian BFGS optimization algorithm.
- Soft line search from the original Nielsen Fortran implementation.
- Trust-region-style control of the line-search input step (`stepmax`).
- Analytic gradients through typed Fortran procedure callbacks.
- Forward and central finite-difference gradients with the upstream step rule
  `abs(x(i))*gradstep(1) + gradstep(2)`.
- Initial inverse-Hessian input in the upstream packed lower-triangle ordering.
- Final inverse Hessian in both packed and full symmetric forms.
- Original convergence codes 1, 2, 3, 4, -2, -4, -5, -6, and -7.
- Gradient-checking computational routine derived from upstream `CHKDFN`.

## Deliberately omitted R-specific code

The `.Call` environment plumbing, R callback bridge, registration code, object
class/name handling, and R printing machinery are not needed in native Fortran.

Upstream `hessian = 1` delegates to the separate R package `numDeriv` rather than
to ucminf's own computational kernel; it is therefore not copied as package
code. The optimizer's BFGS inverse-Hessian estimate is returned directly.

## Modernization

- Free-form Fortran modules and explicit interfaces.
- One package-wide real kind, `dp = real64`, defined once in `ucminf_kinds` and
  re-exported from the public `ucminf` module.
- All real declarations use `real(dp)` and all real literals use `_dp` suffixes.
- No `double precision`, `real*8`, `d0`/`D0` literals, or implicit interfaces.
- Internal packed-symmetric kernels replace external BLAS requirements, keeping
  the package dependency-free for FPM.

## 0.1.1 integration extension

`ucminf_minimize_context` adds a context-bearing callback interface for embedding
UCMINF in other Fortran libraries without global state or nested-procedure
trampolines. The optimization steps are the same as `ucminf_minimize`.
