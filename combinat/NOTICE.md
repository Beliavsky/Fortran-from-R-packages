# NOTICE and provenance

This directory is a modern Fortran translation of computational code from the R package `combinat` version 0.0-8.

Upstream package metadata retained in `upstream/DESCRIPTION` states:

- Package: `combinat`
- Version: `0.0-8`
- Title: `combinatorics utilities`
- Author: Scott Chasalow
- Maintainer: Vince Carey
- License: `GPL-2`
- CRAN publication date: 2012-10-29

The original R files, NAMESPACE, INDEX, and Rd documentation are retained under `upstream/`. Their comments preserve function-specific attribution, including Scott Chasalow, John Wallace, and Alan Zaslavsky.

Algorithm references recorded by the upstream sources include:

- A. Nijenhuis and H. S. Wilf, *Combinatorial Algorithms for Computers and Calculators*, Academic Press, 1978, for combination generation.
- E. M. Reingold, J. Nievergelt, and N. Deo, *Combinatorial Algorithms: Theory and Practice*, Prentice-Hall, 1977, for minimal-change permutations.
- W. Feller, *An Introduction to Probability Theory and Its Applications*, Volume I, 3rd edition, 1968, for generalized binomial coefficients.

The Fortran translation is distributed under GPL-2.0-only, preserving the upstream `GPL-2` licensing designation. The complete GNU GPL version 2 text is included as `LICENSE`.

The translation uses no copied code from external Fortran combinatorics packages and vendors no translated package dependencies. The repository's `rfortran-core` was reviewed before implementation; its shared combinatorial special functions cover integer log-factorial and integer log-choose, but they do not implement the upstream `nCm` operation for real `n`, negative-`n` reflection, or R-style recycling, so no sibling dependency is required here.
