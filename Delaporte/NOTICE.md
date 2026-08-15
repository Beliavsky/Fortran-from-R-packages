# Origin and licensing

This project is a standalone modern-Fortran/FPM adaptation of the computational
code in the R package **Delaporte** version 8.4.3 by Avraham Adler.

Upstream package: https://github.com/aadler/Delaporte

The upstream package declares `BSD_2_clause + file LICENSE`. Its computational
Fortran and C sources carry the BSD 2-Clause terms reproduced in `LICENSE`.
The adapted source files retain copyright/SPDX notices. Copies of the upstream
Fortran computational sources, DESCRIPTION, and CITATION are included under
`upstream/` for provenance.

R interface code, R registration/C wrappers, package metadata machinery, and
thread-control wrappers are not needed by this standalone library. No plotting
code exists in the upstream package.
