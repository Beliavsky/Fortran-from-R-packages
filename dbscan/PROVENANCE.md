# Provenance

## Upstream

- Package: `dbscan`
- Version: 1.2.6
- Date: 2026-08-24
- Upstream URL: https://github.com/mhahsler/dbscan
- Supplied source archive: `dbscan-master.zip`
- Upstream license: GPL (>= 2)

The exact supplied `DESCRIPTION`, `NAMESPACE`, all upstream R sources, and the
package's top-level native C/C++ sources are retained below `upstream/` for
traceability.  The bundled `src/ANN` implementation is intentionally excluded
from the retained source tree to comply with the translation rule against
vendoring dependencies; its license is retained separately in
`licenses/ANN-License.txt`.

## Translation approach

The R/C++ numerical paths were translated into free-form Fortran modules under
`src/`.  A package-specific `dp = real64` kind is defined once in
`dbscan_kinds` and imported throughout.  There is no external linear-algebra or
nearest-neighbor dependency.

The ANN kd-tree backend is replaced by deterministic exact pairwise distance
search.  Cluster labels remain one-based with zero denoting noise, matching the
R package's returned flat labels.  The Fortran result types contain the
computational arrays directly instead of R S3/list attributes.

A check of the target `Fortran-from-R-packages` repository did not identify an
existing top-level `dbscan` translation or a shared dependency providing the
required density-clustering/nearest-neighbor API, so no sibling FPM dependency
is needed for this package.
