# Upstream provenance and license

Source translated from the supplied `mclust-master.zip`.

Upstream package metadata:

- Package: `mclust`
- Version: `6.1.3`
- Date: `2026-07-03`
- Title: *Gaussian Mixture Modelling for Model-Based Clustering,
  Classification, and Density Estimation*
- License: `GPL (>= 2)`

Upstream authors listed in `DESCRIPTION` are Chris Fraley, Adrian E. Raftery,
Luca Scrucca, Thomas Brendan Murphy, and Michael Fop.

Verbatim copies of the supplied `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, and
`inst/CITATION` are in `upstream/`.

## License preservation

The translation is distributed under **GPL-2.0-or-later**, matching the
upstream `GPL (>= 2)` declaration.  `LICENSE` contains the GPL version 2 text;
`LICENSE-GPL-3.txt` is also included because the upstream license permits any
later GPL version.

The original fixed-form numerical files were converted to free-form `.f90` files in `src/legacy/`, with
their upstream comments and provenance intact.  The modern Fortran modules are
a derivative translation and are distributed under the same GPL-2.0-or-later
terms.

BLAS and LAPACK are external numerical libraries and are not copied into this
archive.
