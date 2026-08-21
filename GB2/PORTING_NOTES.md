# Porting notes

## Scope and API

The translation preserves numerical algorithms and exposes them through
array/scalar Fortran interfaces. R graphics, device management, and plotting
wrappers are omitted. R-specific return lists are represented by explicit
arrays and `optimization_result` values.

## Intentional corrections and robustness changes

1. **`main.gb2` / `main2.gb2` RMPG argument typo.** Upstream calls the
   relative median poverty-gap helper with `shape2` as both its third and
   fourth shape arguments. The standalone `rmpg.gb2` definition shows that the
   fourth argument is `shape3`. The Fortran port uses `shape3`.

2. **Moment existence at the boundary.** A GB2 moment requires
   `shape2 + k/shape1 > 0` and `shape3 - k/shape1 > 0`. Equality is divergent.
   The port treats equality as non-existent rather than allowing a gamma pole.

3. **Density at zero.** The port evaluates the mathematical GB2 boundary
   limit: zero when `shape1*shape2 > 1`, finite when it equals one, and
   positive infinity when it is below one.

4. **Quantile ordering.** Scalar `qgb2` is inherently order-independent. The
   compound quantile uses a bracketed Newton/bisection solve against the
   mixture CDF rather than upstream sort/permutation bookkeeping. This avoids
   incorrect reordering for unsorted probability arrays.

5. **Positive-parameter fitting.** Full/profile GB2 fits optimize logarithms of
   positive parameters. This targets the same pseudo-likelihood while avoiding
   invalid negative BFGS iterates.

6. **Compound CDF integration.** The upstream package uses
   `cubature::adaptIntegrate`. This port uses an internal adaptive 15-point
   Gauss-Kronrod integrator with configurable tolerance.

7. **Gini hypergeometric calculation.** The Thomae transformation is preserved,
   but the final convergent real `3F2(1)` is evaluated by an internal series
   implementation instead of `hypergeo::genhypergeo_series`.

8. **Empirical indicators.** `main.emp` no longer depends on R `laeken`; the
   weighted quantile, Gini, mean, ARPR, RMPG, and quintile-share calculations
   are implemented directly. Exact behavior in pathological tied/zero-weight
   samples can differ from `laeken`'s presentation conventions.

9. **Numerical derivatives.** Analytic derivatives are retained where upstream
   supplies them. Indicator Jacobians and a few fitting helper gradients use
   local finite differences, avoiding a direct package-level `numDeriv`
   dependency.

## Dependency licensing

The supplied `cubature-fortran` (`GPL-3.0-or-later`) and
`hypergeo-fortran` (`GPL-2.0-only`) ports cannot both be statically combined
under a common GPL version. They are retained as reference-only snapshots.
The necessary numerical kernels are implemented inside GB2-fortran instead.
See `LICENSES.md`.

## Thread safety

Several optimizers and the compound adaptive-integral callbacks use saved
module context to satisfy procedure interfaces without nonstandard closures.
Those specific operations are not re-entrant across simultaneous threads.
Ordinary density/CDF/moment/indicator evaluations that do not enter these
context-backed routines are otherwise independent.
