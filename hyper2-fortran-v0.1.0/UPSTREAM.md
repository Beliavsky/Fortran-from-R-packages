# Upstream provenance

This project translates the computational code of the R package `hyper2`
version 3.2-3, by Robin K. S. Hankin.

Upstream metadata declares:

- Package: hyper2
- Version: 3.2-3
- Title: The Hyperdirichlet Distribution, Mark 2
- License: GPL (>= 2)

The original package material used for the translation is retained in
`upstream/`, including DESCRIPTION, NAMESPACE, CITATION where present, R source,
and the original Rcpp/C++ sparse-map implementation.

The R package depends on `cubature` and `partitions` for computational tasks.
This release vendors the previously translated pure-Fortran FPM versions of
those packages under `dependencies/`.
