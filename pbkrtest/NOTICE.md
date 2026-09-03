# NOTICE and provenance

This directory is a Fortran translation of computational code from the R
package **pbkrtest**, version 0.5.5.

Upstream package:
- Title: *Parametric Bootstrap, Kenward-Roger and Satterthwaite Based Methods for Test in Mixed Models*
- Authors/copyright holders: Ulrich Halekoh and Søren Højsgaard
- Upstream license: GNU General Public License, version 2 or later
- CRAN publication date recorded by the supplied source: 2025-07-18
- Upstream project URL: <https://people.math.aau.dk/~sorenh/software/pbkrtest/>

The supplied upstream `DESCRIPTION`, `NAMESPACE`, `NEWS`, and `CITATION` files
are retained verbatim in `upstream_metadata/` for attribution and provenance.
The source archive used for this translation was supplied as
`pbkrtest-master.zip`.

## Academic citation

The upstream package asks users to cite:

Ulrich Halekoh and Søren Højsgaard (2014), "A Kenward-Roger Approximation and
Parametric Bootstrap Methods for Tests in Linear Mixed Models -- The R Package
pbkrtest", *Journal of Statistical Software*, 59(9), 1-30.

## Translation scope

The maintained Fortran source translates portable numerical kernels. R object
systems, formulas, model updating/refitting, printing, data-frame/tidy methods,
parallel-cluster orchestration, and R random-number plumbing are deliberately
not copied. The numerical APIs take arrays, scalars, or callbacks directly.

## Reused dependencies

No dependency source is vendored in this directory. The FPM manifest refers to
sibling top-level packages:

- `../rfortran-core` for `dp` and probability distributions;
- `../rfortran-linalg` for SVD, rank, eigensolvers, linear solves, and matrix
  inversion; this shared package owns the pinned `fortran-lapack` dependency;
- `../numDeriv` for callback-based numerical Hessians and Jacobians.

Although the R package depends on `lme4`, its model-object adapter layer is not
translated here. Consequently this numerical package does not add an unused
Fortran `lme4` dependency; callers can combine these kernels with the existing
top-level `lme4` translation when model fitting is required.
