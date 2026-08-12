# Translation notes

## Upstream

- R package: `clue`
- Upstream version: 0.3-68
- Package license field: `GPL-2`
- Main authors: Kurt Hornik; Walter Boehm (contributor)

The complete attached source is retained at `original/clue-master/`.

## Design

The R package uses a flexible S3 framework in which partitions, hierarchies,
memberships and ensembles can be supplied through many third-party R classes.
The Fortran translation replaces this dynamic layer with explicit numerical
representations:

- integer class-id vectors for hard partitions;
- `real(dp)` row-stochastic membership matrices for soft partitions;
- symmetric distance/ultrametric matrices for proximities and hierarchies;
- rank-3 arrays for fixed-size ensembles of membership or hierarchy matrices.

`module clue` reexports the numerical API from the individual implementation
modules.

## Native C translation

`src/assignment.c`, `src/lsap.c`, and the relevant assignment logic are
represented by `clue_lsap.f90`.

The numerical routines from `src/trees.c` are represented by
`clue_trees.f90`.  The triple/quadruple update rules used by iterative
projection and iterative reduction are translated directly at the algorithmic
level rather than replaced by a generic tree package.

## R-level numerical translation

The R-level PAVA, agreement/dissimilarity formulas, medoid models, consensus
updates, validity measures, target ultrametric fits and SUMT logic are
translated to explicit Fortran routines.

For soft Manhattan consensus, the original R routine solves a row-wise linear
program for the L1 fit of a stochastic membership vector.  The Fortran port
uses the supplied `lpSolve-fortran` dependency for the same formulation.
Likewise, exact k-medoids and transportation-based dissimilarities retain their
LP/MILP formulation rather than being replaced with unrelated heuristics.

## Tree-fitting details

The high-level SUMT wrappers use the translated `sumt_optimize` routine and the
translated native penalty functions/gradients.  Unlike the R wrapper's default
multiple random starts, the Fortran convenience functions use the supplied
dissimilarity as a deterministic initial point.  Callers needing multiple
starts can invoke the routines repeatedly with perturbed data/start logic or
use `sumt_optimize` directly.

The R source's penalty-gradient expression adds
`2 * sum(pmin(d, 0))` to every component.  The Fortran high-level SUMT wrappers
preserve that expression for source compatibility.

The target ultrametric fitters accept an `hclust`-style `merge(n-1,2)` matrix:
negative entries identify leaves and positive entries identify earlier merge
rows.  Each topology block is fit by a weighted mean or weighted median and
PAVA enforces nondecreasing merge heights.

## Fortran-specific fixes

R's `&&` and `||` operators short-circuit; Fortran `.and.` and `.or.` are not
required to.  Guarded array accesses were therefore translated to explicit
nested conditions where necessary.  Bounds-checking found and fixed two such
cases during development: the PAVA merge loop and candidate-bin/selection
loops.

Hard partition labels are canonicalized before constructing a membership
matrix.  This is important because R factors/classes need not already be
numbered `1:k`.

## Dependency and licensing

`vendor/lpsolve-fortran` is the user-supplied Fortran translation of lpSolve.
It remains LGPL-2.0-only and has its own metadata, license and original sources.
The top-level clue translation remains GPL-2.0-only.

## Not translated

The numerical library deliberately omits R S3 infrastructure, dynamic method
registration, plotting, data-frame/list presentation, and adapters to external
R clustering packages.  See the README for the extended methods deferred from
this first release.
