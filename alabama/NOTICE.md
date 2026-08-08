# Notices

This is a modern Fortran computational translation of the R package **alabama**
(version 2025.1.0), by Ravi Varadhan with contributions from Gabor Grothendieck.
The original package is licensed under GPL (>= 2). This translation is distributed
under GPL-2.0-or-later.

The bundled `numDeriv-fortran` dependency is the user-supplied Fortran translation
of R package `numDeriv` and retains its GPL-2.0-or-later notices.

The bundled `roptim` dependency is used as a native replacement for the R
`stats::optim` inner solver. It retains its GPL-2.0-or-later notices; its embedded
L-BFGS-B kernel retains the separate BSD-style notice included in that package.
