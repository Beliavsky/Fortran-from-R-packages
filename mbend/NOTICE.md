# NOTICE and provenance

This directory is a Fortran translation of computational code from **mbend
1.3.1**, authored by Mohammad Ali Nilforooshan.

Upstream repository: https://github.com/nilforooshan/mbend

The translated algorithm comes from `R/bend.R`. Upstream documentation credits:

- Jorjani, H. et al. (2003), *A Simple Method for Weighted Bending of Genetic
  (Co)variance Matrices*, Journal of Dairy Science 86:677-679,
  DOI 10.3168/jds.S0022-0302(03)73646-7.
- Schaeffer, L. R. (2014), *Making covariance matrices positive definite*.

## License preservation

The supplied upstream `DESCRIPTION` says `License: GPL-3`, but the supplied
upstream `LICENSE` file is the text of the GNU Lesser General Public License
version 3. This apparent upstream inconsistency has not been silently changed:

- the original `DESCRIPTION` is retained under `upstream/DESCRIPTION`;
- the original `LICENSE` is retained verbatim at both `LICENSE` and
  `upstream/LICENSE`;
- translated Fortran source carries `SPDX-License-Identifier: GPL-3.0-only`,
  following the package metadata's declared license;
- `COPYING` contains the GNU GPL version 3 text for the translated sources.

No upstream copyright, authorship, or citation claim is removed.

## Shared dependency

The translation calls `rfortran-linalg` through the sibling FPM dependency
`../rfortran-linalg` for symmetric eigendecomposition. No copy of that package,
BLAS, or LAPACK is included here.
