# Licensing and provenance

## scip-fortran translation

This translation is distributed under the Apache License 2.0, matching the
attached R package `scip` 1.10.0-3 (`License: Apache License (>= 2)`).
The full Apache-2.0 text is in `LICENSE`.

## Vendored SCIP Optimization Suite

The attached R package vendors SCIP 10.0.2 under `original/scip-master/inst/scip`.
Its `LICENSE` file is retained unchanged there.  SCIP is built as the numerical
backend; the solver source itself has not been rewritten into Fortran.

## Vendored SoPlex

The attached package vendors SoPlex 8.0.2 under
`original/scip-master/inst/soplex`.  Its Apache-2.0 `LICENSE` is retained
unchanged there.

## Original R package

The complete uploaded `scip-master` tree is retained under `original/` for
source-level auditing and attribution.  The package DESCRIPTION identifies
Balasubramanian Narasimhan as author/maintainer and the SCIP Optimization Suite
authors as copyright holders for SCIP, SoPlex, and PaPILO material.
