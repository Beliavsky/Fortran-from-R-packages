# clue-fortran 0.1.0

Modern Fortran/FPM translation of the computational core of R package
`clue` 0.3-68 (Cluster Ensembles).

The original package is by Kurt Hornik with contributions by Walter Boehm and
is distributed under GPL-2.  This translation preserves that license and keeps
the complete attached R/C source tree under `original/clue-master/` for
provenance and algorithm auditing.

## What is translated

### Assignment and partitions

- Linear sum assignment (`solve_lsap`) using a rectangular Hungarian/potential
  algorithm, including maximum-weight assignment.
- Hard and soft membership representations.
- Canonicalization of arbitrary hard class labels.
- Class-id/membership conversion and co-membership matrices.
- Partition lattice meet and join.
- Partition margins, partition coefficient and partition entropy fuzziness.
- Membership construction from cross-dissimilarities.
- Contingency tables.

### Agreements and dissimilarities

Hard-partition agreements include Rand, adjusted Rand, normalized mutual
information, Katz-Powell, Fowlkes-Mallows, Jaccard, purity and prediction
strength.  Membership-based agreements include Euclidean, Manhattan, angle and
diagonal agreement with optimal class matching.

Translated partition dissimilarities include Euclidean, Manhattan,
co-membership, symmetric difference, Rand, GV1, Boehm-A/D/E-style measures,
variation of information, Mallows transportation distance and CSSD.  Mallows
and CSSD use the vendored `lpSolve-fortran` dependency.

Hierarchy/ultrametric matrix comparisons include Euclidean, Manhattan,
cophenetic, gamma, Chebyshev, Lyapunov and spectral dissimilarities, plus the
corresponding principal agreement measures.

### Consensus clustering

- DWH sequential consensus with LSAP class alignment.
- Alternating soft/hard Euclidean consensus.
- Alternating soft/hard Manhattan consensus.
- The soft Manhattan update solves the same row-wise L1 stochastic-vector LP
  used by `clue`, through `lpSolve-fortran`.
- Euclidean hierarchy consensus.
- Ensemble medoid selection from a dissimilarity matrix.

### Medoids and prototype clustering

- Medoid selection.
- Exact `kmedoids` MILP formulation through `lpSolve-fortran`.
- PAM-style k-medoids heuristic.
- Euclidean hard/fuzzy prototype clustering (`pclust_euclidean`).

### PAVA and tree fitting

- Weighted mean and weighted median PAVA.
- Native C tree kernels translated to Fortran:
  - ultrametric deviation and gradient,
  - additive-tree deviation and gradient,
  - iterative projection (IP),
  - iterative reduction (IR).
- Ultrametrification and centroid additive-tree fitting.
- High-level SUMT least-squares ultrametric and additive-tree fits.
- L1 SUMT and IRIP ultrametric fitting.
- Least-squares and L1 ultrametric target fitting for a fixed `hclust`-style
  merge topology.
- Alternating least-squares fit of a sum of ultrametrics.

### Other numerical utilities

- Generic sequential unconstrained minimization (`sumt_optimize`).
- Dissimilarity/variance/deviance accounted-for validity measures.
- Silhouette widths.

## lpSolve dependency

The project vendors the supplied `lpSolve-fortran` 0.1.0 translation under

```
vendor/lpsolve-fortran/
```

and references it as an FPM path dependency.  The dependency remains
LGPL-2.0-only and is kept as a separate package; its license and provenance are
not merged into the `clue` module namespace.

## FPM

```sh
fpm build
fpm test
fpm run --example partition_demo
fpm run --example tree_fit_demo
```

The `fpm.toml` manifest is self-contained because the lpSolve dependency is
vendored by path.

## Minimal example

```fortran
program demo
    use clue
    implicit none
    integer :: a(6), b(6)

    a = [1,1,2,2,3,3]
    b = [2,2,1,1,3,3]

    print *, agreement_adjusted_rand(a,b)
    print *, agreement_nmi(a,b)
end program demo
```

## Validation

Eight regression executables cover assignment, partition metrics, consensus,
LP-backed transport and medoids, PAVA, prototype clustering, native tree
kernels, target fits, SUMT/high-level tree fitting and validity measures.

The release was tested with GNU Fortran 14.2.0 using both optimized compilation
and a bounds-checked warning-as-error build.  See `VALIDATION.md`.

## Translation boundary

The goal is a numerical Fortran library, not an implementation of R's object
system.  Therefore the following are intentionally not reproduced in 0.1.0:

- S3 class/coercion/printing/plotting methods and dynamic method registries;
- adapters to third-party R clustering objects (`kmeans`, `Mclust`, `fanny`,
  `flexmix`, `RWeka`, etc.);
- plotting and dendrogram visualization;
- `cl_bag`/`cl_boot` orchestration around arbitrary R clustering callbacks;
- the most specialized consensus methods (GV3, soft/hard symmetric-difference
  consensus and hierarchy majority consensus);
- hierarchy symmetric-difference and generic Boorman-Olivier integration based
  on R tree class sets;
- the fully generic R `pclust_family` callback/object framework.  The standard
  Euclidean hard/fuzzy numerical family is implemented directly.

These omissions are interface/extended-method scope decisions; the package's
native C assignment and tree-fitting algorithms are translated.
