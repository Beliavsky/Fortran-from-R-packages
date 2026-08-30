# Provenance

## Upstream baseline

This translation targets `spdep` 1.4-2, dated 2026-02-13. Exact package
metadata and the complete upstream contributor list are retained in
`UPSTREAM_DESCRIPTION`; citation metadata is retained in
`UPSTREAM_CITATION.R`.

The implementation was derived by reading the current upstream R algorithms
and the native C kernels in `src/`. It is not a wrapper around R and does not
require an R installation.

## Principal source mappings

| Fortran area | Upstream computational sources used as behavioral references |
| --- | --- |
| Neighbor cardinality, symmetry, components | `src/card.c`, `src/symtest.c`, `src/dfs_ncomp.c`, neighbor utilities in `R/` |
| Distance-band and geodesic neighbors | `src/dnn.c`, `R/dnearneigh.R` |
| K nearest neighbors | `src/knn.c`, `R/knearneigh.R` |
| Gabriel and relative-neighborhood graphs | `src/gabriel.c`, `src/relative.c`, `R/gabrielneigh.R`, `R/relneigh.R`, `R/tri2nb.R` |
| Spatial lags and distance vectors | `src/lagw.c`, `src/nbdists.c`, `R/weights-utils.R`, `R/nbdists.R`, `R/nblag.R` |
| Weight coding and constants | `R/nb2listw.R`, `R/nb2listwdist.R`, `R/autocov.R`, `R/weights-utils.R` |
| Moran statistics | `R/moran.R`, `R/localmoran.R`, `R/moran-bv.R` and native permutation helpers |
| Geary statistics | `R/geary.R`, `R/localC.R`, `src/gearyw.c` |
| Getis-Ord statistics | `R/globalG.R`, `R/localG.R` |
| Join counts | `R/jc.R`, `src/jc.c` |
| Lee association | `R/lee.R`, `R/lee_internal.R` |
| LOSH | `R/LOSH.R` |
| Empirical Bayes and Choynowski | `R/EBI.R`, `R/choynowski.R` |
| Graph utility parity | `R/nboperations.R`, `R/nb2blocknb.R`, `R/rotation.R` |
| MST/SKATER kernels | `R/nbcosts.R`, `R/skater.R`, `src/skater.c` |
| Weighted multivariate spatial delta | `R/delta.R` |

## Fortran-specific design choices

R `nb` and `listw` S3 objects are represented by the `neighbor_list` and
`spatial_weights` derived types. Region indices are one-based, matching both R
and conventional Fortran indexing. Empty neighbor vectors represent isolates;
there is no R-style sentinel value.

The translation uses one package-wide real kind, `dp = real64`, defined once in
`spdep_kinds` and re-exported by the public `spdep` module. All translated
floating-point declarations use `real(dp)` and all real literals use `_dp`
suffixes.

The WGS84 point distance follows the same ellipsoidal approximation and
constants as upstream `src/dnn.c`. KNN is implemented directly in Fortran to
keep this package dependency-free. The Delaunay neighbor routine is a direct
empty-circumcircle implementation; in co-circular degeneracies it can retain
all valid Delaunay edges instead of reproducing `deldir`'s particular tie
triangulation.

The IPFP spatial-delta weights use alternating row/column scaling directly in
Fortran instead of calling R package `mipfp`. Graph-distance weights use an
internal unweighted Floyd-Warshall shortest-path calculation instead of
calling `igraph`. No dependency source has been copied into this package.

The Monte Carlo tests use Fortran's intrinsic pseudorandom generator with an
explicit deterministic seed mapper. A seed therefore gives repeatable results
within a compiler/runtime, but the permutation stream is not expected to be
bit-for-bit identical to R's RNG.
