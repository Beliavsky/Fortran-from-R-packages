# Provenance and licensing

This project is a clean modern-Fortran translation of the computational routines in the CRAN package `cbinom` version 1.6.

Upstream metadata:

- Package: cbinom
- Author: Dan Dalthorp
- Version: 1.6 (2021-04-28)
- License: GPL (>= 2)
- Reference: Andreii Ilienko (2013), "Continuous counterparts of Poisson and binomial distributions and their properties."

This translation is therefore distributed under GPL-2.0-or-later.

The R/Rcpp API mechanics, plotting example, and R-specific vector recycling/NA semantics are not reproduced. The probability algorithms and d/p/q/r functionality are translated to standard Fortran.

The `LICENSE` file contains the GNU GPL version 2 text. The upstream grant permits version 2 or any later version.
