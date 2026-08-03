# Porting notes

## R objects

The R package returns an S3 `hclust` list. The Fortran translation uses the
`hclust_result` derived type. Formula dispatch, labels, calls, and R class
metadata are not meaningful in a standalone Fortran library and are omitted.

## Algorithm selection

The upstream C++ implementation is highly optimized and selects among:

- a minimum-spanning-tree algorithm for single linkage,
- nearest-neighbor chains for reducible linkage methods,
- specialized generic algorithms for centroid and median linkage,
- stored-vector algorithms that avoid materializing all pairwise distances.

The initial Fortran translation uses one deterministic full-matrix
Lance-Williams implementation. This keeps the implementation compact,
self-contained, compiler-portable, and straightforward to validate.

Complexity is therefore:

- memory: `O(n^2)`
- worst-case time: `O(n^3)`

This is a functional translation, not yet a performance-equivalent rewrite of
the specialized C++ kernels.

## Linkage conventions

The formulas are the formulas used by the upstream package:

- single: `min(d(A,C), d(B,C))`
- complete: `max(d(A,C), d(B,C))`
- average: cluster-size weighted average
- McQuitty: equal average
- Ward: Lance-Williams minimum-variance update
- centroid: UPGMC update
- median: WPGMC update

For `ward.D2`, input dissimilarities are squared before updates and merge
heights are square-rooted afterward. For vector Ward clustering, squared
Euclidean distances and the same transformation are used.

For vector centroid and median methods, squared Euclidean distances are updated
and the reported heights are square-rooted, matching the upstream
`hclust.vector` behavior.

## Ties

The engine scans active cluster labels in ascending order and chooses the first
minimum. Valid dendrograms may differ from R, SciPy, or upstream fastcluster
when exact ties permit multiple merge orders. Heights and non-ambiguous merges
remain equivalent.

## Missing values

The distance layer follows the current R-interface formulas in
`src/fastcluster_R.cpp`. It accepts IEEE NaNs as missing coordinates and applies
R's dimension rescaling. The package's legacy pre-R-3.5 Canberra variant is not
exposed because current R uses the modern formula.

## Labels

R row names and `dist` labels are metadata, not numerical computation. The
Fortran result uses integer observation numbers. Applications can keep their
own parallel label array.
