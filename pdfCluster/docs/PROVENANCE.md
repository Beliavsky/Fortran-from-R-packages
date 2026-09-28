# Provenance

## Upstream

Source package: R package `pdfCluster` 1.0-4, published 2022-12-02.

Upstream package metadata and selected source files are retained under
`upstream/`:

- `upstream/DESCRIPTION`
- `upstream/NAMESPACE`
- `upstream/R/pdfCluster.R`
- `upstream/src/kepdf.c`
- `upstream/src/macon.c`
- `upstream/src/nc.c`
- `upstream/inst/CITATION`
- `upstream/inst/CHANGELOG`
- `upstream/NEWS.md`

These files document the R-level algorithms, native numerical kernels, authors,
and release history. They are provenance material only and are not part of the
FPM source tree.

## Translation structure

- `pdfcluster_kinds`: one package-wide `dp = real64` definition.
- `pdfcluster_types`: typed numerical result containers.
- `pdfcluster_stats`: internal numerical helpers.
- `pdfcluster_kde`: `h.norm`, `hprop2f`, and `kepdf` numerical operations.
- `pdfcluster_graph`: valley measure, graph construction, connected components,
  density profiles, and cluster-core extraction.
- `pdfcluster_cluster`: high-level clustering, pair-graph refitting, staged
  classification, and stage extraction.
- `pdfcluster_dbs`: density-based silhouette calculations.
- `pdfcluster_metrics`: adjusted Rand index.
- `pdfcluster`: public facade module.

## External dependency

The translation uses the sibling `geometry` Fortran package solely for
Delaunay triangulation. The expected repository layout is:

    Fortran-from-R-packages/
      geometry/
      pdfCluster/

`pdfCluster/fpm.toml` therefore uses `geometry = { path = "../geometry" }`.
No source from `geometry` is included in this package.

The current sibling `geometry` translation declares `GPL-3.0-or-later`, while
upstream `pdfCluster` declares `GPL-2`. The `pdfCluster` translation preserves
its upstream license declaration; it does not attempt to relicense either
project. Distribution of a combined/link-time work should be reviewed for
license compatibility.
