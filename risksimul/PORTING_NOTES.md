# Porting notes

## Array orientation

Fortran uses assets by simulations for normal inputs, matching the mathematical
layout of upstream `ReturnCopula`. Strata are processed in the same row-major
logical order produced by `as.vector(t(M))` in R.

## Generalized-hyperbolic inverse CDF

The R package constructs `Runuran` PINV generator objects. The Fortran port is
self-contained: it computes a wide GH density grid once per marginal, integrates
that grid, and interpolates the inverse CDF during simulation. The default grid
has 2049 points and can be changed with `gh_grid_size` in `new_portfolio`.

This is deterministic and fast during large simulations, but individual GH
quantiles need not be bit-for-bit identical to `Runuran`.

## Direction optimization

Upstream `Alg3` passes a lower-bound vector whose documented length appears one
shorter than the optimization variable. The Fortran implementation uses a
bounded derivative-free pattern search directly on normalized directions in
the positive orthant. `Alg2` and its rare-event boundary equation are
preserved.

## Orthogonal completion

The upstream routine solves a sequence of small systems using the last matrix
row. The Fortran version uses modified Gram-Schmidt against coordinate vectors.
It guarantees the supplied direction is the first column and checks the final
basis numerically.

## Sufficient statistics

The R implementation stores every weighted observation separately in every
stratum. The Fortran implementation stores counts, sums, squared sums, and
cross-products. This produces the same sample means, unbiased variances, and
covariances while reducing memory from O(total simulations) to O(strata ×
thresholds).

## Allocation corrections

Two upstream edge cases are handled explicitly:

1. `OptAllocHeur` records the previous objective when accepting a better
   candidate. The corrected default records the candidate objective. Set
   `upstream_allocation_compatibility=.true.` to reproduce the literal update.
2. For a positive `mintype` with tail-probability optimization, upstream code
   refers to `s2r`, which only exists in the conditional-excess branch. The
   Fortran implementation selects the requested threshold from the active
   variance objective.

## Minimum stratum size

Like the R package, every generated batch enforces a minimum number of samples
per stratum. Therefore `samples_used` can exceed `samples_requested`, especially
when the requested batch is small relative to the number of strata.

## Conditional excess

The original finite-sample bias-corrected ratio formula is preserved. Both the
unconditional excess numerator and conditional excess are returned, although
the exported R wrappers suppress the unconditional excess in their usual
presentation.

## Random numbers

A deterministic xorshift/Box-Muller/gamma generator is included. Equal seeds
reproduce Fortran results, but streams differ from R.
