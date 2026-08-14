# Upstream and licensing

## rmoo

The supplied source is `rmoo` version 0.3.2, dated 2026-05-03.  Its DESCRIPTION
identifies the license as `GPL (>= 2)` and authors Francisco Benitez and Diego
P. Pinto-Roa.  It describes rmoo as a fork of Luca Scrucca's GA package with
Deb-style nondominated sorting algorithms.

The original `DESCRIPTION`, `NAMESPACE`, and R source directory are retained
under `original/` for provenance and comparison.

## GA dependency

The user supplied `GA-fortran-v0.1.0`, a GPL-2.0-or-later translation of the GA
package.  rmoo only needs its RNG/utilities and basic genetic-operator modules,
which are retained under `src/vendor_ga`.  The original GA DESCRIPTION is
retained under `original/GA/`.

This combined translation is distributed under GPL-2.0-or-later.  `LICENSE`
contains GPL version 2; the upstream "or later" grants permit later GPL versions.
