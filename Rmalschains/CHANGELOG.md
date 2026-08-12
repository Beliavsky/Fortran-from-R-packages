# Changelog

## 0.1.0

* Initial modern Fortran/FPM translation of Rmalschains 0.2-11's R-exposed
  MA-LS-Chains computational path.
* Added SSGA with NAM-3, BLX-alpha, BGA mutation, and worst replacement.
* Added cumulative EA/LS effort scheduling using the original `calculateFrec`
  equation.
* Added persistent per-individual local-search chains.
* Added Solis-Wets, sparse Solis-Wets, persistent simplex, MTS1, MTS2, and
  persistent CMA-ES local search.
* Added target, threshold, initial-population, and package-style accounting.
* Added `actual_nfe` to expose objective calls omitted by upstream counters.
* Added a compatibility switch for the upstream `lsOnly` zero-start/fake-zero
  fitness behavior.
* Corrected the undefined/uninitialized SSW delta clamp to its evident intended
  operation.
* Added five regression executables and two examples.
* Preserved the complete original source tree and license notices.
