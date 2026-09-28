# Notices and provenance

This directory is a modern Fortran translation of computational portions of the R package **ftsa**, version 6.7 (2026-03-31), *Functional Time Series Analysis*.

Upstream authors are Rob Hyndman and Han Lin Shang. The upstream `DESCRIPTION` identifies Han Lin Shang as copyright holder and maintainer and licenses the package under GPL-3. The original `DESCRIPTION`, `NAMESPACE`, and R source files are retained under `upstream/` for provenance and function-by-function comparison.

The Fortran translation is a derivative work distributed under GPL-3.0-only. See `LICENSE`.

The translation uses sibling FPM dependencies `rfortran-core` and `rfortran-linalg` from the Fortran-from-R-packages repository. Dependency source code is not copied into this package. `rfortran-linalg` supplies the checked high-level eigensolver/SVD API used by this translation; this package does not link directly to system BLAS or LAPACK.

R-specific plotting, interactive behavior, S3 object construction, formula parsing, and presentation code are intentionally outside the translation scope. Material numerical compatibility differences are recorded in `fpm.toml`, `README.md`, and `docs/API_COVERAGE.md`.
