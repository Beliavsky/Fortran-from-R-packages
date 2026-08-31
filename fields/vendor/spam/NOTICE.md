# Notices and attribution

This is a translation/port of computational code from **spam 2.11-4**, authored by Reinhard Furrer, Florian Gerber, Roman Flury and contributors. The upstream metadata and citation files are retained under `upstream/`.

The source has multiple upstream licensing/provenance layers. **Do not treat the archive as if every file had one uniform license.** In particular:

- The upstream DESCRIPTION declares `LGPL-2 | BSD_3_clause + file LICENSE`.
- `upstream/0LICENSE` separately states that the R code and documentation are under the GPL, that most other Fortran routines are from Youcef Saad's SparseKit under LGPL, and records the provenance of the modified sparse-Cholesky routines.
- The Ng-Peyton/PCx Cholesky notice reproduced verbatim in `upstream/0LICENSE` grants no-charge use/reproduction/derivative works/redistribution subject to retention and documentation conditions, but separately says commercial organizations or incorporation into a product for sale require contacting the named Argonne representative. That notice is retained unchanged and should be reviewed for any intended commercial distribution/use.
- Bundled ARPACK sources retain the Rice University copyright identified by the package `LICENSE`; the corresponding BSD-3-Clause text is in `LICENSES/BSD-3-Clause-ARPACK.txt`.
- The user-supplied `r_mod.f90` was supplied as MIT-licensed code. The exact original is retained as `upstream/r_mod-original.f90`; its build copy only changes free-form wrapping/continuation layout.

Key numerical contributors credited by the upstream package include Youcef Saad (SparseKit); Esmond G. Ng, Barry W. Peyton, Joseph W.H. Liu and Alan D. George (sparse Cholesky/ordering); and Rich Lehoucq, Kristi Maschhoff, Danny Sorensen and Chao Yang (ARPACK).

Recommended spam citation retained upstream includes:

Reinhard Furrer and Stephan R. Sain (2010), *spam: A Sparse Matrix R Package with Emphasis on MCMC Methods for Gaussian Markov Random Fields*, Journal of Statistical Software 36(10), 1-25, DOI 10.18637/jss.v036.i10.

This NOTICE summarizes provenance for convenience; the verbatim upstream files control where they differ from this summary.
