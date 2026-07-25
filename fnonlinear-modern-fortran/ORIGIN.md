# Origin and licensing

This project was translated from the attached `fNonlinear` source package:

- Package: `fNonlinear`
- Version: `4052.83`
- Package date: 2026-02-21
- Declared license: `GPL (>= 2)`

The original metadata and R/C sources used for the translation are retained in
`original/` for provenance.

Several numerical components in the R package are themselves derived from
older GPL-compatible work:

- BDS code by Blake LeBaron, adapted for R by Adrian Trapletti
- Nonlinearity tests adapted from the `tseries` package
- Nonlinear dynamics kernels derived from `tseriesChaos`/TISEAN work by
  Antonio Fabio Di Narzo and collaborators
- Fixed-step RK4 code historically distributed with `odesolve`

All new Fortran source files carry `GPL-2.0-or-later` SPDX identifiers and an
explicit GNU GPL version 2-or-later notice. `LICENSE` contains the GNU GPL
version 2 text.
