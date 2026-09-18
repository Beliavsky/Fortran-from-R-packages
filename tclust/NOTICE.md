# NOTICE and provenance

This directory is a modern free-form Fortran translation of computational functionality from the R package **tclust**, upstream version 2.2-1.

Upstream project: <https://github.com/valentint/tclust>

Upstream package authors and contributors identified by `DESCRIPTION` include Valentin Todorov, Luis Angel García Escudero, Agustín Mayo Iscar, Javier Crespo Guerrero, and Heinrich Fritz. The original package is distributed under GPL-3. The included `LICENSE` contains the GNU General Public License version 3 text. Copyright, authorship, references, and package provenance remain with the upstream authors as applicable.

The `upstream/` directory retains selected computational R/C++ sources plus upstream metadata solely to make provenance and algorithm correspondence auditable. Plotting and presentation-only R sources are not copied into this translation.

## Numerical implementation notes

The Fortran implementation is dependency-free. It does not vendor Rcpp, RcppArmadillo, Armadillo, BLAS, LAPACK, `r.f90`, `r_mod.f90`, or another translated R package. Symmetric eigendecompositions use a package-local Jacobi implementation, which is adequate for the small covariance matrices used by the translated algorithms and avoids a system BLAS/LAPACK requirement.

The translation follows the upstream algorithms for trimmed Gaussian clustering, trimmed k-means, robust linear grouping, eigenvalue and determinant/shape restrictions, discriminant diagnostics, information criteria, simulations, and agreement indices. Exact R RNG-stream matching is not attempted. The GPCM covariance-pattern restriction machinery in `src/GPCM.cpp` is retained for provenance but is not yet translated into the public Fortran fitter.

No code from an external dependency has been copied into the Fortran source tree.
