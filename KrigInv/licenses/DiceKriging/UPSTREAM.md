# Upstream provenance

Translated source package:

- Package: DiceKriging
- Version: 1.6.1
- Date: 2025-10-20
- Title: Kriging Methods for Computer Experiments
- Authors listed by the supplied `DESCRIPTION`: Olivier Roustant, David
  Ginsbourger, Yves Deville, with Clément Chevalier and Yann Richet as
  contributors.
- Upstream license declaration: `GPL-2 | GPL-3`

The original `DESCRIPTION`, `NAMESPACE`, and `CHANGELOG` from the supplied
archive are retained under `upstream/`.

The translated covariance formulas were checked against the package's native C
sources `src/CovFuns.c` and `src/Scaling.c`, while likelihood, LOO, kriging,
SCAD, CV, and update behavior were translated from the corresponding R source
files.

The numerical regression targets in the permanent Fortran tests are derived
from the upstream package's own `tests/testthat/test-km.R` and related test
expectations. No R executable is required to run the Fortran tests.
