# Translation notes

## Upstream

This project translates the computational portion of CRAN package `cubature`
2.1.4-1.  That package wraps two independent native libraries:

1. Steven G. Johnson's `cubature` 1.0.4 (`hcubature` and `pcubature`).
2. Thomas Hahn's Cuba 4.2.2 (`Cuhre`, `Divonne`, `Suave`, and `Vegas`).

The original R wrappers and native sources used for this translation are kept
under `upstream/` for provenance.

## What was translated

### Johnson cubature library

`hcubature` is a direct algorithmic translation of the important numerical
parts of `hcubature.c`:

- 15-point Gauss-Kronrod / embedded 7-point Gauss rule in one dimension;
- embedded Genz-Malik degree-7/degree-5 fully symmetric rule in dimensions
  greater than one;
- local error estimation;
- split-dimension selection from fourth-difference indicators;
- global refinement by repeatedly bisecting the region with the largest
  estimated error;
- individual, paired, L1, L2, and Linf convergence criteria.

The Fortran implementation uses array storage rather than the C binary heap.
That changes region-management complexity and exact evaluation ordering, but
not the cubature rule or refinement criterion.

`pcubature` retains the upstream p-adaptive tensor-product Clenshaw-Curtis
method.  The original C code has an elaborate nested-value cache and can
increase polynomial degree anisotropically.  This Fortran v0.1 implementation
recomputes each nested tensor grid and raises the degree isotropically.  It is
therefore target-equivalent but can use more evaluations, especially in four
or more dimensions.

### Cuba methods

The R package's complete public computational method surface is exposed, but
v0.1 does not reproduce every Cuba C optimization line-for-line:

- `cuhre` uses the translated deterministic h-adaptive embedded cubature
  engine.  It targets the same deterministic integral and error contract, but
  does not reproduce Cuba's key-selected Cuhre rule tables internally.
- `divonne` uses adaptive region partitioning with randomized low-discrepancy
  estimates.  `key1` controls the baseline regional sample size.  Cuba's
  extrema search, `xGiven`, and peak-finder subprotocol are not reproduced in
  v0.1.
- `suave` combines adaptive region subdivision with randomized
  low-discrepancy sampling and independent regional error estimates.  It
  preserves the method's adaptive stratification role but not Cuba's exact
  fluctuation/grid formulas.
- `vegas` has a native adaptive separable importance grid, independent
  sampling within bins, inverse-variance combination over iterations, and
  iterative grid rebinning.  State-file persistence and named reusable Cuba
  grids are omitted.

Thus the Johnson `hcubature` path is the closest source-level translation;
the Cuba paths preserve the public numerical targets and distinct algorithm
families but are not bit-for-bit ports of Cuba 4.2.2.  The original Cuba C
sources are included under `upstream/Cuba/` to make these differences easy to
audit and to support a future exact-backend port.

## Vector interface

Fortran scalar and vector callback interfaces are both provided.  In v0.1 the
vector entry points preserve numerical semantics but may invoke a vector
callback with a single point in some algorithm paths; the upstream C/R code
can batch hundreds of points for performance.  This affects throughput, not
the target integral.

## Infinite limits

`cubintegrate` implements the same tangent mapping used in the R wrappers:

`x = tan(t)`, with Jacobian `sec(t)^2`, component by component.

This supports finite, semi-infinite, and doubly infinite boxes.

## Result-code and probability differences

The upstream Cuba routines report a method-specific failure count and a
chi-square reliability probability for each component.  In v0.1 the Fortran
result uses a uniform return-code convention (`0` success, `1` evaluation or
region limit, `2` bad argument, `3` internal failure), and `prob(:)` is set to
`1` as a placeholder.  Integral and error estimates are the validated
computational outputs; exact Cuba chi-square reliability bookkeeping is a
remaining fidelity item.

For infinite-bound `cubintegrate` calls, a module procedure pointer is used
only during the synchronous tangent-transformed integration.  Consequently,
that particular wrapper is not reentrant/thread-safe in v0.1; the finite-bound
method entry points are unaffected.

## Omitted R-only infrastructure

The following are intentionally not translated:

- Rcpp registration and `.Call` glue;
- R list/class construction and argument matching;
- verbose R console output;
- package vignette/benchmark infrastructure;
- fork-worker controls tied to the Cuba C process model;
- Cuba state-file serialization and R callback metadata (`cuba_weight`,
  `cuba_iter`, `cuba_phase`).

No plotting code is part of the numerical translation.

## Validation

The packaged tests cover:

- the upstream `prod(cos(x))` 2-D result;
- the Wang-Landau 1-D reference value `1.63564436296`;
- exact polynomial integration;
- p-adaptive Clenshaw-Curtis integration;
- the two-component Cuba phase/shift example;
- all four Cuba-family entry points;
- scalar/vector callback interfaces;
- infinite-bound Gaussian integration.

The second component of the upstream phase/shift test is independently
checked against high-accuracy numerical integration as
`0.3078074096213368...`; the R test file stores the rounded value `0.3078155`.
