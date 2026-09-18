# dbscan — modern Fortran translation

This directory translates the computational core of the R package **dbscan
1.2.6** to modern free-form Fortran with FPM.  The public module is
`dbscan_api`.

The translation covers density-based clustering, exact nearest-neighbor and
fixed-radius search, OPTICS, shared-nearest-neighbor methods, LOF, a substantial
HDBSCAN implementation, DBCV, and related graph/structural helpers.  Plotting,
tidying, printing, and other R-specific presentation code is intentionally not
translated.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The package has no external Fortran dependency.  In particular, it does **not**
vendor or link the upstream ANN library.  Neighbor searches use a deterministic
exact brute-force backend.  This preserves exact Euclidean neighbor semantics
and makes the package self-contained on Windows, at the cost of O(n^2) search
work instead of ANN's kd-tree acceleration.

## Public API

Import the public surface with:

```fortran
use dbscan_api
```

The main translated procedures include `knn`, `frnn`, `dbscan`, `optics`,
`extract_dbscan`, `extract_xi`, `snn`, `snnclust`, `jpclust`, `lof`,
`pointdensity`, `hdbscan`, `coredist`, `mrdist`, `dbcv`, `ncluster`, and
`nnoise`.  Derived result types retain numerical results without reproducing R
S3 classes or attributes.

`hdbscan` supports the upstream `cluster_selection_epsilon` HDBSCAN(e)
criterion and returns GLOSH scores from the upstream condensed-hierarchy
noise-exit/death-distance calculation. `extract_fosc` supports unsupervised
extraction plus symmetric full-matrix pair constraints, mixed `alpha`
objectives, and unstable-branch pruning. R adjacency-list/dist-vector constraint
forms and R object augmentation are intentionally outside the Fortran API.

## Translation coverage

Package status: **substantial**.  Coverage is **25 of 25 (100%)** exported
computational R functions.  This fraction measures functions with a meaningful
Fortran mapping; it does **not** mean complete R-interface or behavioral parity.
The intentionally partial entries are called out below and in `fpm.toml`.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `adjacencylist` | `dbscan_nn, dbscan_api` | `adjacencylist_knn, adjacencylist_frnn` | substantial |
| `as.dendrogram` | `dbscan_conversion` | `optics_to_dendrogram` | partial |
| `as.reachability` | `dbscan_conversion` | `dendrogram_to_reachability` | partial |
| `comps` | `dbscan_api` | `comps_knn, comps_frnn` | substantial |
| `coredist` | `dbscan_hdbscan` | `coredist` | substantial |
| `dbcv` | `dbscan_dbcv` | `dbcv` | substantial |
| `dbscan` | `dbscan_cluster` | `dbscan, dbscan_from_frnn` | substantial |
| `extractDBSCAN` | `dbscan_optics` | `extract_dbscan` | complete |
| `extractFOSC` | `dbscan_hdbscan` | `extract_fosc` | substantial |
| `extractXi` | `dbscan_optics` | `extract_xi` | substantial |
| `frNN` | `dbscan_nn` | `frnn, frnn_query, frnn_from_dist` | substantial |
| `glosh` | `dbscan_hdbscan` | `glosh` | substantial |
| `hdbscan` | `dbscan_hdbscan` | `hdbscan` | substantial |
| `is.corepoint` | `dbscan_cluster` | `is_corepoint` | substantial |
| `jpclust` | `dbscan_snn` | `jpclust` | substantial |
| `kNN` | `dbscan_nn` | `knn, knn_query, knn_from_dist` | substantial |
| `kNNdist` | `dbscan_nn` | `knn_dist` | substantial |
| `lof` | `dbscan_outlier` | `lof` | substantial |
| `mrdist` | `dbscan_hdbscan` | `mrdist` | substantial |
| `ncluster` | `dbscan_api` | `ncluster` | complete |
| `nnoise` | `dbscan_api` | `nnoise` | complete |
| `optics` | `dbscan_optics` | `optics` | substantial |
| `pointdensity` | `dbscan_outlier` | `pointdensity` | substantial |
| `sNN` | `dbscan_snn` | `snn` | substantial |
| `sNNclust` | `dbscan_snn` | `snnclust` | substantial |

Notable remaining gaps are R `dist`/S3 dispatch and object construction,
FOSC adjacency-list/dist-vector constraint forms and automatic constraint
repair, exact general R dendrogram/reachability object conversion, and
ANN-specific approximate/kd-tree search controls. See `API_COVERAGE.md` for
details.

## Validation

`test/test_dbscan.f90` contains deterministic regression tests for all major
algorithm families, including exact GLOSH reference scores, constrained/mixed
FOSC extraction, and the upstream package's 118-point HDBSCAN(e) regression
case (`cluster_selection_epsilon = 1`), which must yield five clusters and no
noise. `validation/` also contains an independent comparison with scikit-learn
for DBSCAN partitioning, LOF values, HDBSCAN partitioning, and OPTICS
reachability/core distances. HDBSCAN membership probabilities are validated
against the R package's documented core-distance formula because scikit-learn
uses a different definition.

See `VALIDATION.md` for the exact compiler flags and environment limitations.

## Licensing and provenance

The upstream R package is GPL (>= 2).  `COPYING` contains GPL version 2; the
translation is distributed under GPL-2.0-or-later subject to all retained
upstream notices.  The upstream package bundles ANN, whose separate license is
preserved at `licenses/ANN-License.txt`; ANN source is **not** included or used
by this translation.  See `NOTICE.md` and `PROVENANCE.md`.
