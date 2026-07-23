# Fortran from R packages

This repository collects experimental modern Fortran translations and ports of
computational code from several R packages. Each subdirectory is an independent
Fortran Package Manager (fpm) project with its own documentation, tests,
provenance record, and license.

These projects are unofficial and are not endorsed by the original package
authors, CRAN, or the R Foundation. They are experimental and incompletely
validated; independently verify numerical results before relying on them.

## Packages and licenses

| Directory | Original R package | License for the Fortran project |
| --- | --- | --- |
| `fGarch-modern-fortran/` | fGarch 4052.93 | GPL-2.0-or-later |
| `longmemo-modern-fortran/` | longmemo 1.1-4 | GPL-2.0-or-later |
| `rugarch-modern-fortran/` | rugarch 1.5-6 | GPL-3.0-only |
| `timsac-modern-fortran/` | timsac 1.3.8-6 | GPL-2.0-or-later |
| `tseries-modern-fortran/` | tseries 0.10-62 | GPL-2.0-only OR GPL-3.0-only |

The repository is an aggregate of separately licensed projects; there is no
additional repository-wide license. Consult the license/copying, notice,
origin, and README files within each package before redistributing or modifying
it.

## Attribution and provenance

Each project identifies the original package, version, authors or copyright
holders, and the nature of the Fortran work. Original package metadata is
retained under package `reference/` directories where available. The modern
Fortran files also carry SPDX license identifiers and dated modification
notices.
