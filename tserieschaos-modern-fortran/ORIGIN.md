# Origin and licensing

This project is a clean modern Fortran translation of the computational routines in:

- Package: `tseriesChaos`
- Version: 0.1-13.1
- Title: Analysis of Nonlinear Time Series
- Original author and maintainer: Antonio Fabio Di Narzo
- Package date: 2013-04-29
- CRAN publication date recorded in the supplied source: 2019-01-07
- Original license declaration: `GPL-2`

The supplied `DESCRIPTION` and `NAMESPACE` files are retained under `original/` for provenance.

The translation uses the SPDX expression `GPL-2.0-only`, includes the complete GNU GPL version 2 text in `LICENSE`, and places a GPL-2.0-only notice in every Fortran source file.

The original package notes that its nonlinear-time-series routines were largely inspired by the TISEAN project. The algorithm references from the original documentation remain applicable, including Hegger, Kantz, and Schreiber (1999), Kennel, Brown, and Abarbanel (1992), Eckmann, Kamphorst, and Ruelle (1987), and Provenzale et al. (1992).

The supplied C sources `search.c`, `find_knearests.c`, and `false.nearest.c` define the original box-assisted neighbor search. Version 0.2.0 translates that computational path into a dynamic modern Fortran box index while retaining the direct implementation as a reference mode.
