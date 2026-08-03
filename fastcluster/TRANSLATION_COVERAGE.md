# Translation coverage

Upstream package: `fastcluster` 1.3.0.

## Exported R routines

| Upstream routine | Fortran representation | Status |
|---|---|---|
| `hclust` | `hclust`, `hclust_matrix`, `hclust_condensed` | translated |
| `hclust.vector` | `hclust_vector` | translated |

## Linkage methods

| Method | Status |
|---|---|
| single | translated |
| complete | translated |
| average/UPGMA | translated |
| mcquitty/weighted/WPGMA | translated |
| ward.D | translated |
| ward.D2 | translated |
| centroid/UPGMC | translated |
| median/WPGMC | translated |

## Vector metrics

| Metric | Status |
|---|---|
| Euclidean | translated |
| maximum | translated |
| Manhattan | translated |
| Canberra | translated |
| binary | translated |
| Minkowski | translated |

## Internal computational behavior

- Lance-Williams distance updates: translated
- R-compatible merge numbering: translated
- Dendrogram leaf ordering: translated
- Initial membership sizes: translated
- Pairwise missing-coordinate handling: translated
- Infinite dissimilarities: supported
- Input-preservation behavior: inherent through `intent(in)`

## Deliberately not translated

- R `.Call` registration and protection code
- S3 class construction and print/plot integration
- R/Python installation infrastructure
- LaTeX vignette build system
- Specialized C++ MST, nearest-neighbor-chain, heap, and stored-vector
  performance kernels

The omitted optimized kernels affect asymptotic performance, not the available
linkage methods or numerical update formulas. The original implementation is
retained in `original/fastcluster-master/src/`.
