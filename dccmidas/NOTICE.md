# NOTICE

This directory is a modern Fortran translation of computational portions of the R package **dccmidas 0.1.3**.

The upstream package is authored and maintained by **Vincenzo Candila** and is distributed under the GNU General Public License version 3. The original package metadata, selected R sources, native C++ wrappers, citation file, and references are retained under `upstream/` for attribution and numerical provenance.

The translation does not vendor the source of the package's R dependencies. It is designed to reuse sibling packages in `Fortran-from-R-packages`, in particular `rfortran-core`, `rfortran-linalg`, `roll`, `rumidas`, `rugarch`, and `maxLik`.

The upstream native helpers `Det` and `Inv` use RcppArmadillo/Armadillo. The Fortran translation does not copy or compile Rcpp, RcppArmadillo, or Armadillo; it maps those operations to the shared `rfortran-linalg` API.

No BLAS, LAPACK, R compatibility layer, or translated dependency source is copied into this package.

See `PROVENANCE.md` for the source-to-source mapping and material compatibility differences.
