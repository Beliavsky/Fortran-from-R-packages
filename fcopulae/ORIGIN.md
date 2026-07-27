# Origin and licensing

This project was translated from the attached source archive:

- Package: `fCopulae`
- Original version: `4052.86`
- Original date: `2026-02-21`
- Original license metadata: `GPL (>= 2)`

The original `DESCRIPTION` and `ChangeLog` are retained as `ORIGINAL_DESCRIPTION` and `ORIGINAL_CHANGELOG`.

The modern Fortran translation is distributed under **GPL-2.0-or-later**. Every Fortran source, application, example, and test file contains an SPDX identifier and an explicit GPL version 2-or-later notice. `LICENSE` contains the full GNU GPL version 2 text.

The package originally depends on R packages such as `fBasics` and `fMultivar`. Required probability, integration, optimization, RNG, and linear-algebra pieces are included locally so the Fortran library does not call R.
