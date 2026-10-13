# Upstream provenance and notices

This package translates CRAN **CreditRisk 0.1.7**, published 2024-04-19,
by Alessandro Cimarelli (author/maintainer) and Nicolò Manca (author).

- Package page: https://cran.r-project.org/package=CreditRisk
- Release source: https://cran.r-project.org/src/contrib/CreditRisk_0.1.7.tar.gz
- Archive: https://cran.r-project.org/src/contrib/Archive/CreditRisk/
- Target repository: https://github.com/Beliavsky/Fortran-from-R-packages
- Original license: MIT + file LICENSE; copyright year 2017, copyright
  holder Alessandro Cimarelli & Nicolò Manca, preserved byte-for-byte in
  `upstream/LICENSE`.

`upstream/DESCRIPTION`, `upstream/NAMESPACE` and `upstream/R/*.R` preserve the
original metadata, computational source, documentation comments, author
attribution and references. `upstream/SHA256SUMS` identifies those exact
source files. They are reference material only, not compiled or vendored
third-party dependency implementations. Original plotting examples and
R-specific checks remain in the reference copies, but are not translated.
The binary `cdsdata` example dataset is not redistributed.

The original DESCRIPTION cites Damiano Brigo, Massimo Morini and Andrea
Pallavicini (2013), *Counterparty Credit Risk, Collateral and Funding*,
DOI 10.1002/9781118818589; Print ISBN 9780470748466, Online ISBN
9781118818589, © 2013 John Wiley & Sons Ltd. Other upstream references include
David Lando (2004), *Credit Risk Modeling*; E. G. Haug, *The Complete Guide to
Option Pricing Formulas*; and C. Hastings (1955), *Approximations for Digital
Computers*. Original detailed references and simulation-book authors remain
in the preserved R source comments. No text from those books is newly copied.

Fortran changes include explicit typed interfaces, pure model computations,
validation/status results, exact limiting cases, stable log-tail multiplication,
stable exponential integrals, reproducible local simulation, and linked
L-BFGS-B calibration with diagnostics. See README/API for material source
formula quirks and intentional differences. In particular the incorrect
BlackCox reporting of initial L/K is corrected, while Merton's survival-rate
term and CDS's accrued-premium subtraction are preserved and explicitly tested.

## Reused dependencies

Before implementation, existing top-level sibling APIs were inspected.
- `rfortran-core`: `r_kinds` (dp, i64, pi), `r_distributions` (normal CDF and
  log-CDF), `r_transforms` (stable expm1). Manifest license MIT; authored by
  Vivek Rao and contributors.
- `lbfgsb3`: objective-only `lbfgsb_minimize_fd`, control/result types and
  L-BFGS-B 3.0 implementation. Manifest license GPL-2.0-only. Its own NOTICE
  and LICENSES preserve optimizer/authors' notices. Follow its GPL terms when
  distributing combined linked binaries; this port's MIT license does not
  override them.

Both use sibling FPM path dependencies. No dependency code, BLAS, LAPACK,
ARPACK, r.f90 or r_mod.f90 is copied here. No system numerical-library link
flags or direct LAPACK dependency are needed. Shared dependency sources are
unchanged. The port includes no R runtime or translated stats dependency.
