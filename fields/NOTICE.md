# Notices and attribution

This project is a Fortran translation/port of computational code from **fields 17.3**, authored by Douglas Nychka, Reinhard Furrer, John Paige, Stephan Sain, Florian Gerber, Matthew Iverson, Rider Johnson, and contributors.

The upstream package metadata declares `License: GPL (>= 2)`. The original metadata, citation, R sources, C sources, and fixed-form Fortran numerical source used for this translation are retained under `upstream/`. The translated `fields` source is distributed under GPL-2.0-or-later.

The upstream package asks users to cite:

Douglas Nychka, Reinhard Furrer, John Paige and Stephan Sain, *fields: Tools for spatial data*, University Corporation for Atmospheric Research. The upstream package documentation gives DOI 10.5065/D6W957CT and project URL `https://github.com/dnychka/fieldsRPackage`.

## Bundled spam dependency

The user supplied `spam-fortran-v0.1.0`, which is bundled under `vendor/spam` as an FPM path dependency. Its licensing is not uniform. Its `NOTICE.md`, retained upstream notices, and `LICENSES/` directory control the provenance of that dependency. In particular, its own notice records separate terms/provenance for SparseKit, Ng-Peyton/PCx Cholesky code, ARPACK, and `r_mod`.

A small local compatibility correction was made in `vendor/spam/src/spam_cholesky.f90`: when `pivot='none'`, the permutation and inverse permutation are explicitly initialized to the identity before calling the inherited Ng-Peyton `cholstepwise` routine. The supplied wrapper previously passed an uninitialized permutation on this branch, causing a segmentation fault. This changes wrapper initialization only; the inherited sparse-Cholesky numerical routine is otherwise unchanged.

The upstream `fields/LICENSE.note` itself warns that spam contains Fortran routines with distinct licensing/provenance issues. Review the spam notices for any intended redistribution or commercial use.

This NOTICE is a convenience summary. Verbatim upstream files govern where they differ.
