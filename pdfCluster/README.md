# pdfCluster

Modern free-form Fortran translation of the computational parts of the R package
`pdfCluster` 1.0-4, which implements density-based clustering using multivariate
kernel density estimation.

The translation is intended for the root of
`Beliavsky/Fortran-from-R-packages` as the top-level directory `pdfCluster`.
It uses the sibling `geometry` translation for Delaunay triangulation:

```toml
[dependencies]
geometry = { path = "../geometry" }
```

No Qhull, Rcpp, BLAS, LAPACK, or system numerical library is copied or linked by
this package.

## Build

From the `pdfCluster` directory, with the sibling `geometry` directory present:

```text
fpm build
fpm test
fpm run --example basic
```

The public convenience module is `pdfcluster`.

## Numerical scope

The translation includes:

- Gaussian and Student-t7 product-kernel density estimation;
- fixed and adaptive bandwidths and the normal-reference bandwidth rule;
- unidimensional, Delaunay, and pairwise-valley neighborhood graphs;
- valley-area calculation and graph connected components;
- density-profile construction and con2tree-style cluster-core extraction;
- staged density-ratio classification;
- density-based silhouette diagnostics; and
- adjusted Rand index calculation.

The Delaunay path calls the public `delaunayn_points` API of the sibling
`geometry` package. The pairwise graph uses a modern Fortran translation of the
upstream valley-measure algorithm.

## Translation coverage

Package status: **substantial**.

**8 of 8 (100.0%)** exported computational R functions are mapped. This fraction
measures mapped computational functions, not complete R/S4 interface
compatibility. Plotting, printing, summaries, class declarations, and other
presentation-only methods are excluded from the denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `kepdf` | `pdfcluster_kde` | `kepdf_fit`, `kepdf_estimate` | substantial |
| `h.norm` | `pdfcluster_kde` | `h_norm` | complete |
| `hprop2f` | `pdfcluster_kde` | `hprop2f` | complete |
| `pdfCluster` | `pdfcluster_cluster`, `pdfcluster_graph` | `pdfcluster_fit`, `pdfcluster_regraph_pairs`, `density_profile`, `core_labels_from_profile` | substantial |
| `pdfClassification` | `pdfcluster_cluster` | `pdf_classification` | substantial |
| `groups` | `pdfcluster_cluster` | `groups_at_stage` | substantial |
| `dbs` | `pdfcluster_dbs` | `dbs_from_labels`, `dbs_from_pdfcluster` | substantial |
| `adj.rand.index` | `pdfcluster_metrics` | `adjusted_rand_index` | complete |

Important interface differences include the omission of R S4 object construction,
R data-frame coercion, plotting and dendrogram creation, arbitrary R bandwidth
function lookup, and Qhull option strings. Integer zero is used where R objects
would store `NA` for observations not yet assigned to a cluster. The translated
adaptive classification path avoids cumulative in-place multiplication of the
core bandwidths across stages.

See `docs/API_COVERAGE.md` for details and `fpm.toml` for the machine-readable
mapping.

## License and provenance

The upstream R package declares `License: GPL-2`; the translation preserves
that declaration and includes the GNU GPL version 2 text in `LICENSE`. Original
R/C source excerpts used for provenance are retained under `upstream/` and are
not compiled into the Fortran package.

The sibling Fortran `geometry` package currently declares `GPL-3.0-or-later`.
Because this package preserves the upstream `GPL-2` declaration, downstream
redistributors of a combined/link-time work should review the licenses and their
compatibility rather than assuming the two declarations can be relicensed or
combined without restriction. See `NOTICE.md` and `docs/PROVENANCE.md`.
