# Translation coverage

## Translated computational algorithms

- Variable-shape Nelder-Mead reflection/expansion/outside contraction/inside contraction/shrink.
- Spendley-Hext-Himsworth fixed-shape reflection/reflection-next/shrink method.
- Box constrained complex method, including reflection line search, bound projection,
  and scaling toward a feasible reference point for nonlinear positive inequalities.
- Initial simplexes: given, axes, Spendley regular simplex, Pfeffer/MATLAB-style,
  and random-within-bounds.
- Simplex sorting, center, centroid excluding vertices, size, function spread,
  variance, standard deviation, shrink, and simplex-gradient estimate.
- MATLAB-like `fminsearch` and bounded `fminbnd` front ends.
- Direct `fmin_gridsearch` over `npts^n` grid points.
- Iteration/function-evaluation termination, simplex-size termination,
  size-plus-function-spread termination, Box repeated-match termination,
  variance/std termination, and Kelley stagnation detection.
- O'Neill and Kelley restart detection and automatic restart loops.
- History storage and typed output/user-termination callbacks.

## Dependency handling

The R implementation delegates its bookkeeping and simplex data structure to the
R packages `optimbase` and `optimsimplex`. This Fortran project translates the
subset of their *computational roles needed by neldermead* directly into native
Fortran modules instead of recreating their R object systems and string-key
get/set APIs.

## Deliberate API differences

- R list/S3 objects and `neldermead.set`/`neldermead.get` string-key plumbing are
  replaced by `nm_options`, `nm_result`, and `nm_simplex` derived types.
- Objective and constraint functions use explicit Fortran procedure interfaces.
  Nonlinear inequalities follow the upstream convention `c(x) >= 0`.
- R plotting, shell logging/formatting, `OutputFcn`/`PlotFcns` R-list marshalling,
  and S3 `print`/`summary` methods are not translated.
- The R-only `method='mine'` extension hook is replaced by directly writing a
  Fortran optimizer or using the provided typed callbacks; it is not a numerical
  algorithm belonging to the package.
- `optimbase.gridsearch` is reimplemented locally according to the documented
  `fmin.gridsearch` range rule.
- The exact internal `optimsimplex(method='oriented')` restart construction is
  unavailable in the supplied package; the Fortran oriented restart contracts
  the previous simplex by 1/2 about its best vertex, preserving orientation.

## Constraint fidelity

Box feasibility uses the upstream convention that every nonlinear inequality is
nonnegative. Candidate points are projected to bound constraints and, when needed,
moved toward a feasible reference point with `box_ineq_scaling` until feasible or
`gui_alpha_min` is reached.

For random-bound Box initialization the Fortran port first accepts independently
sampled feasible vertices. If a sampled vertex is infeasible it tries bounded
rejection sampling before falling back to the upstream-style geometric scaling
toward the feasible starting/reference point. If geometric scaling reaches
`gui_alpha_min` without obtaining a numerically feasible point, the reference point
itself is retained. This is a robustness adjustment for nearly active nonlinear
constraints: it prevents a stochastic initial complex from collapsing almost
entirely onto `x0`, while preserving feasibility and the Box search equations.

The reflection line search likewise enforces bounds and nonlinear feasibility on
the reflected point before repeated contraction toward the centroid. This ordering
is slightly more defensive than the R implementation around floating-point
constraint boundaries and is documented here rather than hidden.
