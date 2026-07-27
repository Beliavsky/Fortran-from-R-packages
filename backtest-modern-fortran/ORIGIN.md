# Origin and licensing

This project is a clean modern Fortran translation of the computational
routines in the CRAN package:

- Package: `backtest`
- Version: 0.3-4
- Date: 2015-09-17
- Authors: Jeff Enos, David Kane, and contributors listed in the original
  `DESCRIPTION`
- Original license declaration: `GPL (>= 2)`

The translation preserves that license as `GPL-2.0-or-later`.

The original source archive supplied for translation was
`backtest-master.zip`. No R source or data files are required to build the
Fortran project. The algorithms were translated from the R source files under
`R/`, especially `backtest.function.R`, `backtest.compute.R`,
`calc.turnover.R`, `overlaps.compute.R`, `categorize.R`, `bucketize.R`, and
the numerical parts of `backtest.R`.

The original authors do not endorse this translation. Numerical equivalence is
claimed only for the definitions and test cases documented in `VALIDATION.md`.
