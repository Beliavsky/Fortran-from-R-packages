# Upstream provenance

- Package: `ranger`
- Version: 0.18.0
- Upstream package date: 2026-01-15
- CRAN publication recorded in DESCRIPTION: 2026-01-16
- Uploaded source archive used for this translation: `ranger-master.zip`
- Source archive SHA-256: `221ec3e7eac4871a1af718a16f1a3e48eea55adae8141464676ca28b22730e51`

The exact per-file SHA-256 hashes for upstream `R/` and `src/` files are recorded in `upstream/SHA256SUMS.txt`. The archive hash is repeated in `upstream/ARCHIVE_SHA256.txt`.

The translation was made from the attached source tree rather than from generated CRAN binaries. Plotting/presentation/R-object infrastructure was intentionally excluded. Statistical computations implemented in R and C++ were ported to typed modern Fortran APIs.

Dependency review for the target `Fortran-from-R-packages` repository found no existing top-level `ranger` translation. The translation reuses the repository's focused `rfortran-core` package for `r_kinds::dp` and does not copy that package into this directory.
