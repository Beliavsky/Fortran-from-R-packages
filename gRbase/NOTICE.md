# Notices and provenance

## Upstream package

This directory is a clean Fortran translation of computational algorithms from
**gRbase 2.0.3**, an R package by Søren Højsgaard. The attached upstream
snapshot identifies the package as GPL (>= 2), i.e. GPL-2.0-or-later.

Upstream metadata is retained verbatim in:

- `upstream/DESCRIPTION`
- `upstream/CITATION`

The primary upstream citation is Claus Dethlefsen and Søren Højsgaard (2005),
"A Common Platform for Graphical Models in R: The gRbase Package", Journal of
Statistical Software 14(17), 1-12.

## Combination algorithm attribution

Upstream `src/combn-bristol/combnPrimC.c` states that the combination algorithm
is by R. F. Warren-Smith (private communication), with P. T. Wallace, Starlink,
dated 25 August 1999, and:

> Copyright (C) 1999 Rutherford Appleton Laboratory

That upstream code is licensed under GPL version 2 or, at the recipient's
option, any later version. This translation retains the algorithmic provenance;
the original C source is not copied into the Fortran package.

## Minimal triangulation and maximal-prime decomposition

Upstream `R/igraph_mintriang.R` and `R/igraph_mpd.R` credit Clive Bowsher, with
modifications by Søren Højsgaard. The algorithms refer to the recursive
thinning and aggregate-cliques methods described by Kristian G. Olesen and
Anders L. Madsen (2002), "Maximal Prime Subgraph Decomposition of Bayesian
Networks", IEEE Transactions on Systems, Man, and Cybernetics, Part B, 32(1).
This attribution is retained for `minimal_triangulation` and
`maximal_prime_decomposition`.

## Dependency provenance

This translation reuses existing sibling packages rather than copying their
sources:

- `rfortran-core` for the common `dp` real kind;
- `rfortran-linalg` for SPD linear algebra and its pinned `fortran-lapack`
  dependency;
- `igraph` / FPM package `igraph-fortran` for maximal-clique and connected-
  component algorithms.

No Rcpp, RcppEigen, RcppArmadillo, Eigen, Armadillo, Matrix, BLAS, LAPACK or
igraph dependency source is vendored here.

## Translation notes

R S3/S4 dispatch, Rcpp `.Call` registration, R object coercion/dimnames,
printing, plotting, package datasets, and other R-specific interface machinery
are intentionally outside the Fortran translation. Fortran APIs use explicit
integer labels and arrays instead.
