# Translation notes

## Source package

- R package: matchingMarkets 1.0-5
- License declaration: GPL (>= 2)
- Original package copyright notice: Thilo Klein, 2014-present

The original package also contains code derived from the Stable Matching Suite
(Fahiem Bacchus) and MiniSat-era support code (Niklas Sorensson and others).
Those original files are retained unchanged under `original/`; the new
Fortran algorithms do not copy their C++ implementation text.

## Main design choices

1. Preference matrices are numeric integer matrices. Missing/unacceptable
   choices use `0` rather than R `NA`.
2. Results use allocatable derived types instead of R lists/data frames.
3. The R/Java constraint implementations of `hri` and `sri` are represented by
   native exact enumeration. This preserves the mathematical solution set on
   tractable instances, not the original SAT/constraint-engine performance.
4. `hri3` follows the round/batch semantics in `src/eadam.cpp`, including the
   efficiency-adjustment outer loop and consent-gated interrupting pairs.
5. `hri2` uses an exact bounded-search couples solver in v0.1.0. This is a
   correctness/reference path for small and medium instances, not a
   performance replacement for the original Roth-Peranson/Bacchus matcher.
6. `plp` formulates the same binary partitioning LP as the R function and
   solves it through the supplied `lpSolve-fortran` translation.
7. KHB is exposed at matrix level with an integer `z_index` instead of R
   column-name/formula manipulation.
8. `stabit` and `stabit2` are deferred rather than replaced with a different
   selection model.

## Fortran portability

All procedure interfaces are explicit. Index-guard logic does not depend on
short-circuit evaluation of `.and.` or `.or.`. Source lines fit the standard
132-column free-form limit.
