# Translation coverage

## Translated

- Rowan Subplex outer iteration (`subplx.f`).
- Nelder-Mead subspace optimizer (`simplx.f`).
- Subspace partitioning (`partx.f`).
- Step adaptation/reorientation (`setstp.f`).
- Initial simplex construction, centroid updates, reflection/expansion/
  contraction/shrink operations, simplex ordering, stable distance, and
  decreasing-key sorting.
- R wrapper semantics for scalar/vector `parscale`, convergence codes, and
  finite-difference Hessian calculation.

## Modernization

- Fortran 77 COMMON state is represented by a private derived type passed
  explicitly between routines.
- BLAS calls used only for copying, scaling, AXPY, and absolute sums are
  replaced by intrinsic array operations; no external BLAS is required.
- The objective callback has an explicit modern Fortran procedure interface.
  Every invocation occurs in a module-level routine with a non-optional typed
  dummy procedure, avoiding host-associated implicit-interface diagnostics on
  older gfortran/FPM combinations.
- Exact coincidence tests from the original code are expressed without real
  `==`/`/=` syntax so strict `-Werror` builds remain warning-clean.

## Omitted R-only code

- `.Call`/SEXP registration and R environment handling.
- Parameter names and `...` argument forwarding; callers use ordinary Fortran
  closures/module state or explicit data structures instead.

There is no plotting code in the upstream package.

## Scope note

The historical continuation/user-option `mode` documented in Rowan's original
routine is not exposed by the current R package; Aaron King's upstream source
hardcodes first-call/default-option behavior. The Fortran translation matches
that current package behavior rather than restoring the old continuation API.
