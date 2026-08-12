# Translation notes

## Full matching

The Fortran `fmatch_core` reproduces the network construction in upstream
`R/fmatch.R`:

1. each eligible treatment-control pair is a unit-capacity arc carrying its
   discrepancy as cost;
2. treatment, control, End, and Sink nodes receive the same capacity/supply
   bookkeeping used by the R implementation;
3. an integral min-cost-flow problem is solved;
4. selected treatment-control arcs are converted into connected matched-set
   labels.

`fullmatch` also reproduces the orientation rule used by the R implementation.
When the permitted control:treatment ratio is predominantly below one, or a
negative omission fraction requests omission from the treatment side, the
problem is transposed and the reciprocal constraints are solved on the flipped
network.

The R package discretizes floating-point discrepancies before sending them to
RELAX-IV in some paths.  This Fortran implementation solves directly with
real(dp) arc costs, so it does not need the R tolerance/discretization layer.
The resulting optimization target is the original discrepancy sum.

## Minimum-cost flow

The internal solver is a successive shortest augmenting-path algorithm on an
integral residual network.  Bellman-Ford is used for residual shortest paths so
reverse arcs with negative cost are handled correctly.  The full-matching
network starts with nonnegative forward costs, matching the upstream distance
requirements.

## Mahalanobis calculations

The ordinary Mahalanobis routine reproduces the pooled within-group covariance
used by upstream `compute_mahalanobis()`:

`S = ((nt-1) S_t + (nc-1) S_c) / (nt+nc-2)`.

If the covariance is singular, a symmetric generalized inverse is formed from
a Jacobi eigendecomposition with a numerical rank threshold.

Rank Mahalanobis follows upstream `src/smahal.cc`: columns are ranked with
average ranks for ties, covariance is computed over ranks, tied marginal rank
variances are rescaled to the untied sample variance `n(n+1)/12`, and the
resulting covariance is generalized-inverted.

## Sparse distances

The R `InfinitySparseMatrix` is represented by `type(sparse_distance)` for
compact storage or by `type(distance_spec)` with an explicit logical eligibility
mask.  This avoids depending on R's S4 slots and makes forbidden pairings
unambiguous even when IEEE infinities are undesirable.

## R-specific functionality

Formula parsing and model fitting are intentionally outside this library.  For
example, the R `glm` method ultimately forms a one-dimensional fitted score and
then computes standardized absolute treatment-control differences.  In
Fortran, compute the model elsewhere and pass the resulting score vector to
`score_distance()`; use `standardization_scale()` when the pooled scaling is
wanted.
