# NOTICE

This project contains a modern Fortran translation of computational algorithms from the R package `ecp` 3.1.6.

Copyright and authorship remain with the respective upstream authors and contributors. Upstream `DESCRIPTION` identifies Nicholas A. James, Wenyu Zhang, and David S. Matteson as authors and declares `GPL (>= 2)`.

The original package implements native routines with Rcpp/C++. Those sources are retained under `upstream/src/` only for provenance and review. The maintained implementation under `src/` is free-form Fortran and does not link to or require Rcpp.

See `docs/PROVENANCE.md` and the per-function mappings in `fpm.toml`.
