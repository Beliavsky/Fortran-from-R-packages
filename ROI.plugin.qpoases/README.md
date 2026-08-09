# ROI.plugin.qpoases-fortran

Modern Fortran/FPM translation of the computational portion of
`ROI.plugin.qpoases` 1.0-3 and its bundled qpOASES 3.2 solver.

The package solves convex quadratic programs

\[
  \min_x \frac12 x^T H x + g^T x
\]

subject to

\[
  lb \le x \le ub,\qquad lb_A \le A x \le ub_A .
\]

The native Fortran implementation uses a feasible primal working-set method.
Equalities are permanently active, blocking bounds/constraints enter the
working set, and constraints with negative KKT multipliers leave it.  A
one-variable auxiliary Phase-I QP is used to obtain a feasible starting point.
Persistent models retain the previous solution and active-set IDs for hotstarts.

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example roi_example
fpm run --example hotstart_example
```

## Main modules

* `qpoases` - solver types, model API, QProblem/QProblemB/SQProblem-like entry
  points, hotstarts, diagnostics, options and status constants.
* `roi_qpoases` - array-level analogue of the ROI plugin conversion layer.
* `qpoases_active_set` - feasible working-set and Phase-I algorithms.
* `qpoases_linalg` - standalone dense linear algebra.

No R, Rcpp, BLAS or LAPACK dependency is required.

## Hessian type constants

`hst_zero`, `hst_identity`, `hst_posdef`, `hst_posdef_nullspace`,
`hst_semidef`, `hst_indef`, and `hst_unknown` retain qpOASES' integer values
0 through 6.

## Licensing

The R plugin declares GPL-3.  Bundled qpOASES 3.2 source is
LGPL-2.1-or-later.  Source files in this translation carry the license
appropriate to their provenance.  Because the public plugin/API layer is
GPL-3-derived, distribute the combined package under GPL-3-compatible terms.
See `NOTICE.md`, `UPSTREAM_PROVENANCE.md`, and `LICENSES/`.

## Scope

This release translates the computational surface exercised by the R package:
general and simply bounded QPs, persistent models, hotstarts, bounds/general
constraints, dual/primal/objective getters, active-set diagnostics, Hessian
types, options storage, ROI conversion, infeasible detection and the important
zero-Hessian LP unbounded case.

It does not claim a line-for-line translation of every advanced qpOASES 3.2
internal.  See `TRANSLATION_COVERAGE.md` for exact differences.
