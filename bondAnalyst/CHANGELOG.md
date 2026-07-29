# Changelog

## 1.0.1-fortran.1 - 2026-07-28

- Ported all 55 exported `bondAnalyst` computational routines to modern
  Fortran.
- Added an FPM project layout, application, example, and two automated test
  programs.
- Added status-code and IEEE-NaN error handling.
- Replaced R polynomial-root filtering with bracketed nonnegative yield and
  spread solvers.
- Added `effectiveannualratezcb` and `conventionalpercentpricechange` while
  preserving the original source formulas under their original names.
- Preserved GPL-3 licensing, copyright attribution, original source, package
  metadata, and references.
