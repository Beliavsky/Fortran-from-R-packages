# Fortran from R packages

This repository collects experimental modern Fortran translations and ports of
computational code from 23 R packages. Each subdirectory is an independent
Fortran Package Manager (fpm) project with its own documentation, tests,
provenance record, and license.

These projects are unofficial and are not endorsed by the original package
authors, CRAN, or the R Foundation. They are experimental and incompletely
validated; independently verify numerical results before relying on them.

## Packages and licenses

| Directory | Original R package | License for the Fortran project |
| --- | --- | --- |
| `bayesgarch-modern-fortran/` | bayesGARCH 2.1.10 | GPL-2.0-or-later |
| `betategarch-modern-fortran/` | betategarch 3.4 | GPL-2.0-only |
| `deoptimr-modern-fortran/` | DEoptimR 1.2-0 | GPL-2.0-or-later |
| `fbasics-modern-fortran/` | fBasics 4052.98 | GPL-2.0-or-later |
| `fbonds-modern-fortran/` | fBonds 3042.78 | GPL-2.0-or-later |
| `fcopulae-modern-fortran/` | fCopulae 4052.86 | GPL-2.0-or-later |
| `fextremes-modern-fortran/` | fExtremes 4032.84 | GPL-2.0-or-later |
| `fGarch-modern-fortran/` | fGarch 4052.93 | GPL-2.0-or-later |
| `fmultivar-modern-fortran/` | fMultivar 4031.84 | GPL-2.0-or-later |
| `fnonlinear-modern-fortran/` | fNonlinear 4052.83 | GPL-2.0-or-later |
| `garchx-modern-fortran/` | garchx 1.7 | GPL-2.0-or-later |
| `gogarch-modern-fortran/` | gogarch 0.7-6 | GPL-2.0-or-later |
| `lgarch-modern-fortran/` | lgarch 0.7 | GPL-2.0-only |
| `longmemo-modern-fortran/` | longmemo 1.1-4 | GPL-2.0-or-later |
| `msgarch-modern-fortran/` | MSGARCH 2.51 | GPL-2.0-or-later |
| `rmgarch-modern-fortran/` | rmgarch 1.4-2 | GPL-3.0-only |
| `robustbase-modern-fortran/` | robustbase 0.99-7 | GPL-2.0-or-later |
| `rugarch-modern-fortran/` | rugarch 1.5-6 | GPL-3.0-only |
| `stochfactor-modern-fortran/` | stochvol 3.2.9 and factorstochvol 1.1.2 | GPL-2.0-or-later |
| `timsac-modern-fortran/` | timsac 1.3.8-6 | GPL-2.0-or-later |
| `tsdyn-modern-fortran/` | tsDyn 11.0.5.2 | GPL-2.0-or-later |
| `tseries-modern-fortran/` | tseries 0.10-62 | GPL-2.0-only OR GPL-3.0-only |
| `tserieschaos-modern-fortran/` | tseriesChaos 0.1-13.1 | GPL-2.0-only |

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
