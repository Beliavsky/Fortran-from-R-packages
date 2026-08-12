# Translation notes

## Upstream

- Package: matchingR
- Version: 2.0.0
- Date in DESCRIPTION: 2025-09-22
- Authors: Jan Tilly and Nick Janetos
- License: GPL (>= 2)

## Native-code mapping

`src/galeshapley.cpp` maps to `src/matchingr_galeshapley.f90`. The Fortran
translation preserves the package's strict reviewer comparison (`>`), which
means ties do not poach a currently held proposal.

`src/roommate.cpp` maps to `src/matchingr_roommate.f90`. Phase 1 proposal
holding, table reduction, rotation discovery/elimination, and the odd-market
dummy used by the R wrapper are retained. Ragged C++ deques are represented by
bounded Fortran preference-list objects.

`src/toptradingcycle.cpp` maps to `src/matchingr_ttc.f90`. The implementation
uses the equivalent cycle-elimination formulation: each active agent points to
its highest-ranked active owner, every directed cycle is finalized, and the
process repeats. TTC's core outcome is unique under strict complete
preferences, so the result is independent of which cycle is discovered first.

`src/utils.cpp` maps to `src/matchingr_utils.f90`.

## College admissions

matchingR implements capacities by expanding a college into repeated slots and
then applying one-to-one Gale-Shapley. The Fortran implementation follows the
same construction for both student-optimal and college-optimal modes, including
heterogeneous `slots(:)`.

## Index conversion

matchingR's C++ core uses IDs `0..n-1`; the R interface converts most results to
`1..n`. Fortran public results are always 1-based with zero as an unmatched
sentinel. Preference-input functions accept either convention.

## Roommate checker bugs

The original checker has two observable scan issues, and the R wrapper also
mixes 0-based preferences with 1-based returned matchings. These bugs are not
used by the Irving search itself. The Fortran checker is corrected by default;
`legacy_cpp_bug=.true.` preserves the two C++ source-level scan quirks in a
1-based translation for audit/comparison purposes.

## Omitted code

- Rcpp generated registration glue
- R list/matrix/NA marshalling
- warnings and printing
- documentation/vignette machinery

The complete upstream tree is retained unchanged under `original/`.
