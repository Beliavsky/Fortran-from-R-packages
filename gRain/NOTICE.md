# NOTICE and provenance

## Upstream package

This work translates computational code from the R package **gRain**, version
**1.4.6**, published on CRAN on 2026-03-02.

Upstream author and maintainer:

- Søren Højsgaard

Upstream description: probability propagation in Bayesian networks (graphical
independence networks).

The original package declares `License: GPL (>= 2)`. This translated package is
therefore distributed under **GPL-2.0-or-later**. The GNU GPL version 2 text is
included as `LICENSE` and `COPYING`; the "or later" permission comes from the
upstream license declaration retained in `upstream/DESCRIPTION`.

The upstream `DESCRIPTION` and `inst/CITATION` files are retained verbatim under
`upstream/`. The recommended citation identifies:

Søren Højsgaard (2012), "Graphical Independence Networks with the gRain Package
for R", Journal of Statistical Software 46(10), 1-26,
doi:10.18637/jss.v046.i10.

## Translation provenance

The source snapshot supplied for this translation was the user-provided archive
`gRain-master.zip`. Its compiled layer contains Rcpp/C++ implementations in
`src/propagate.cpp` and a small sparse-matrix helper in `src/toBeRemoved.cpp`.
The numerical functionality of gRain is broader than those C++ files, so the
translation also follows the package's R implementations for network
compilation, evidence handling, queries, simulation, CPT construction, and
categorical-data estimation.

The Lauritzen-Spiegelhalter propagation implementation in
`src/grain_propagation.f90` is a native Fortran translation of the collect and
distribute logic in upstream `src/propagate.cpp`, expressed through gRbase
probability-table primitives.

## Reused sibling packages

Before implementation, the public `Fortran-from-R-packages` repository was
checked for an existing gRain translation and compatible shared modules. No
existing gRain translation was found. The translation intentionally reuses the
sibling `gRbase` package and `rfortran-core` rather than copying their source.
The gRbase sibling in turn owns its graph and linear-algebra dependencies.

No source from gRbase, igraph, rfortran-core, rfortran-linalg, BLAS, LAPACK,
ARPACK, Rcpp, Eigen, or Armadillo is vendored in this directory.

## Interface adaptations

R character node names, factor levels, formulas, S3/S4 object dispatch, and
R-specific attributes are replaced by explicit integer node/state identifiers
and typed Fortran structures. Probability-table values preserve R/Fortran
column-major ordering. CPTs store the child variable first, followed by its
parents.

The Fortran `simulate_network` routine does not reproduce R's RNG stream. It
uses a portable Park-Miller generator; `simulate_from_uniforms` is provided for
caller-controlled random streams. This is an interface/runtime adaptation, not
a change to the conditional sampling probabilities.

Hugin NET parsing/writing, plotting, printing/summary formatting, formula
parsing, broom/predict integration, and other R-facing convenience layers are
not translated because they are interface functionality rather than portable
numerical kernels.

No affiliation with or endorsement by the upstream authors is implied.
