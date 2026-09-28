# NOTICE

## Upstream package

This project is a modern Fortran translation of computational routines from the
R package `pdfCluster` version 1.0-4.

Upstream authors listed in `DESCRIPTION`:

- Giovanna Menardi — author and maintainer
- Adelchi Azzalini — author
- Tiziana Rosolin — contributor through version 0.1-13

The upstream package declares `License: GPL-2`. The GNU General Public License
version 2 text is included as `LICENSE`. The retained files under `upstream/`
are included for attribution and provenance and are not compiled by FPM.

The original native implementation includes C routines for kernel density
estimation, valley measures, and connected-component traversal. The maintained
Fortran sources in `src/` are translations/reimplementations of the numerical
operations and do not compile or vendor the upstream C sources.

## geometry dependency

The R package imports the R package `geometry`. This Fortran translation uses
the sibling Fortran translation through:

    geometry = { path = "../geometry" }

Only the public Delaunay API is used. No geometry source is copied into this
package.

At the time this translation was prepared, the sibling Fortran `geometry`
package declared `GPL-3.0-or-later`, whereas upstream `pdfCluster` declares
`GPL-2`. This project does not purport to relicense either code base. Anyone
distributing a combined/link-time work should independently review the license
compatibility and obtain any permissions that may be required.

## Material compatibility differences

- R S4 objects, plotting, summaries, printing, and dendrogram presentation are
  omitted.
- Delaunay construction uses the sibling pure-Fortran `geometry` implementation
  rather than the R package's Qhull-backed call and does not expose Qhull option
  strings.
- R `NA` cluster assignments are represented by integer zero.
- Arbitrary R functions supplied through `h.funct` are replaced by typed
  numerical entry points.
- The translated adaptive `hcores` classification path avoids cumulative
  in-place bandwidth multiplication across stages.
