# Algorithm notes

## caRamel generation rules

`new_xval` retains the five generation mechanisms used by the R implementation:

1. interpolation in objective-space Delaunay simplexes touching the Pareto front;
2. extrapolation along oriented improvement edges;
3. multivariate normal generation using the reference population covariance structure;
4. parameter/block recombination from Pareto-archive parents;
5. the periodic "fireworks" independent-variable perturbations around individual-objective maxima and the maximin point.

The R implementation's rank standardization constant `a = 3/8`, minimization-to-maximization transformation, covariance amplification `sqrt(2)`, extrapolation boost, exponential step length, and bounds clipping are retained.

## Native n-dimensional Delaunay triangulation

Upstream `newXval` calls `geometry::delaunayn`, which in turn relies on external geometry/Qhull machinery. The Fortran port is self-contained and uses an incremental Bowyer-Watson-style n-dimensional triangulation:

- construct an enclosing super-simplex;
- insert objective points one at a time;
- remove simplexes whose circumspheres contain the new point;
- identify the boundary facets of the resulting cavity;
- retriangulate the cavity with the new point;
- discard all final simplexes touching a super-simplex vertex.

Circumcenters are obtained by pivoted Gaussian elimination. A very small deterministic perturbation is applied internally to resolve exact co-spherical/degenerate ties; returned simplex indices always refer to the original points. The super-simplex is deliberately far from the data (`1000*d*radius`) because randomized differential testing showed that a merely enclosing but less distant simplex could omit occasional convex-hull cells after the artificial vertices were removed.

## Port corrections

Several corrections were made where a literal translation would preserve accidental R behavior rather than the intended algorithm.

### Duplicate objective rows in `newXval`

Upstream calls `unique(obj)` before triangulation but subsequently indexes the original `param`, `crit`, and `Fo` arrays with simplex indices. When a duplicate occurs before a later unique row, these indices can refer to the wrong original individual. The Fortran port maintains an explicit first-occurrence map from each unique objective row back to its original population row.

### Population downsize selection

The R `downsize()` code contains scalar/logical handling around `which(D == 0)` that is fragile when several candidates or extrema are present. The port implements the evident intent: among minimum-front-rank candidates in a box, preferentially retain a candidate that is an absolute maximizer of an objective; otherwise select randomly.

When a front-rank boundary must be sampled to enforce `popsize`, the Fortran port samples without replacement. This prevents duplicate population indices, which weighted sampling with replacement could otherwise introduce.

### Additional objective outputs

The R main loop attempts to carry additional values returned after the optimization objectives, but later resets a matrix dimension using only `nvar+nobj`, which is inconsistent when additional columns exist. The Fortran `nout` option keeps these columns consistently in the archive and final population while optimizing only the first `nobj` values.

### Sensitivity call count

The R sensitivity loop evaluates `nfront*nvar` perturbed points but increments `nrun` only once by `nfront`. The Fortran result reports the actual number of objective callback invocations.

## Pareto duplicates

The upstream Pareto implementation treats a duplicate as dominated by the first retained equal point. The Fortran implementation preserves this convention using deterministic lexicographic ordering with original row order as the final tie-breaker.

## Random numbers

The statistical rules are preserved, but numerical random streams are not expected to match R because the port uses the Fortran intrinsic RNG plus a Box-Muller normal transform. `seed_random(seed)` provides reproducible runs for a given compiler/runtime.
