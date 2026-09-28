# NOTICE and provenance

This package is a modern free-form Fortran translation of the computational
routines in **Rlinsolve 0.3.3**, authored and maintained upstream by Kisung You.
The supplied upstream archive identifies the package as MIT licensed and was
packaged on 2025-09-21 for CRAN publication on 2025-09-22.

The original R package uses RcppArmadillo and the R `Matrix` package. This
translation does not copy Rcpp, Armadillo, Matrix, BLAS, or LAPACK sources and
does not link to system BLAS/LAPACK. The iterative algorithms are implemented
in package-local free-form Fortran and use a small partial-pivoting LU routine
for user-supplied preconditioner solves.

For auditability, the relevant upstream R and C++ computational sources,
DESCRIPTION, NAMESPACE, original R-package LICENSE metadata, and REFERENCES.bib
are retained under `upstream/`. They are provenance material and are not built
by FPM.

Primary upstream references retained in `upstream/inst/REFERENCES.bib` include
Saad, Demmel, Hestenes and Stiefel, van der Vorst, and the Templates for the
Solution of Linear Systems literature cited by Rlinsolve.
