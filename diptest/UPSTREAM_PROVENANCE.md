# Upstream provenance and dependency review

- Package: **diptest 0.77-2**, dated 2025-08-19, CRAN publication 2025-08-20.
- Canonical CRAN page: <https://CRAN.R-project.org/package=diptest>
- Source archive: <https://cran.r-project.org/src/contrib/diptest_0.77-2.tar.gz>
- Source archive SHA-256:
  `3babf6964be9865a3d872c95c8c8ed49953240e7db5717718893a110d3ebcb3f`
- CRAN mirror source revision used for initial inspection:
  `f7806b8f5cffb79efcd191bbdd4282b4ae7b68ec`.
- Original table `inst/extraData/qDiptab.rds` SHA-256:
  `90079b2d7168d5d796a06ae0c51b137d048ee3162df5ac3944ce8e14a03c1705`.
- Translation date: 2026-10-12.

The retained C kernel and RDS table were byte-compared with the official CRAN
archive. `upstream/` retains the original package DESCRIPTION and NAMESPACE,
`R/dip.R`, `R/dipTest.R`, `src/dip.c`, the table, the statfaculty dataset, and
relevant manual pages. These are reference and attribution material. FPM
compiles only the maintained free-form source in `src/`, not the original C.

`tools/generate_table.R` reads the original binary64 table and emits 17-digit
round-trip decimal constants. It retains all 546 values in their original
column-major order, the 21 sample-size nodes, and the 26 probability nodes.
There is no replacement asymptotic approximation or regenerated simulation table.

## Repository dependency check

The target repository was inspected before implementing dependencies at revision
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.

- Reused `rfortran-core` through `path = "../rfortran-core"`.
- `r_kinds: dp` is the single real kind for every maintained Fortran source.
- `r_sorting: r_sort` provides compatible pure sorting of a copied input.
- No compatible dip statistic implementation was identified in the inspected
  target package/shared modules; the package-specific geometry is translated.
- No linear algebra is required, so `rfortran-linalg`, BLAS, LAPACK, ARPACK and
  `rfortran-compat` are not dependencies.
- The intrinsic RNG is explicitly documented as distinct from R's generator;
  deterministic uniform matrices avoid a generator dependency.

No dependency source, binaries, generated build files or dependency caches are
part of the deliverable. The package must be installed as a top-level directory
beside the existing core package, not built in isolation without its dependency.
