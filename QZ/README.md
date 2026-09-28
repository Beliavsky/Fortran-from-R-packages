# QZ - modern Fortran translation

This directory translates the computational API of the R package **QZ 0.2-4** to modern free-form Fortran. The upstream package provides standard and generalized eigendecompositions, Schur/QZ decompositions, and Schur/QZ reordering on real and complex matrices.

The translation calls the modern explicit LAPACK interfaces from the `fortran-lapack` package rather than copying the fixed-form LAPACK wrappers shipped in the R source. The configured dependency is the Beliavsky fork of Federico Perini's `fortran-lapack`; see `NOTICE.md`.

## Build

The package is intended to live as a top-level directory in `Fortran-from-R-packages` and uses the FPM-compatible LAPACK fork directly:

```toml
[dependencies]
fortran-lapack = { git = "https://github.com/Beliavsky/fortran-lapack.git" }
```

Build and test with:

```text
fpm build
fpm test
fpm run --example basic
```

No separately installed BLAS or LAPACK library is required by the package manifest.

## Public Fortran API

`qz_api` re-exports the package API and `dp` kind. The main entry points are:

- low-level wrappers: `qz_dgees`, `qz_dgeev`, `qz_dgges`, `qz_dggev`, `qz_dtgsen`, `qz_dtrsen`, and their complex `z*` counterparts;
- high-level generics: `qz_decompose`, `qz_geigen`, `geigen`, and `ordqz`;
- `fda_geigen` for the augmented generalized-eigenanalysis used by `fda::geigen`;
- `h_transpose` for the Hermitian transpose.

Results are returned in typed derived types rather than R lists/S3 objects. LAPACK `INFO` values and translation-specific dimension errors are retained in each result's `info` field.

## Translation coverage

Package status: **substantial**.

**18 of 18 (100.0%)** exported computational R functions have meaningful Fortran mappings. This percentage measures mapped computational functions, not complete R-interface compatibility. Printing methods, `.onLoad`/`.onUnload`, S3 formatting, R dynamic-library plumbing, and R object construction are intentionally outside the coverage basis.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `qz.dgees` | `qz_lapack` | `qz_dgees` | complete |
| `qz.dgeev` | `qz_lapack` | `qz_dgeev` | complete |
| `qz.dgges` | `qz_lapack` | `qz_dgges` | complete |
| `qz.dggev` | `qz_lapack` | `qz_dggev` | complete |
| `qz.dtgsen` | `qz_lapack` | `qz_dtgsen` | complete |
| `qz.dtrsen` | `qz_lapack` | `qz_dtrsen` | complete |
| `qz.zgees` | `qz_lapack` | `qz_zgees` | complete |
| `qz.zgeev` | `qz_lapack` | `qz_zgeev` | complete |
| `qz.zgges` | `qz_lapack` | `qz_zgges` | complete |
| `qz.zggev` | `qz_lapack` | `qz_zggev` | complete |
| `qz.ztgsen` | `qz_lapack` | `qz_ztgsen` | complete |
| `qz.ztrsen` | `qz_lapack` | `qz_ztrsen` | complete |
| `fda.geigen` | `qz_highlevel` | `fda_geigen` | substantial |
| `ordqz` | `qz_highlevel` | `ordqz` | substantial |
| `qz` | `qz_highlevel` | `qz_decompose` | substantial |
| `qz.geigen` | `qz_highlevel` | `qz_geigen` | substantial |
| `geigen` | `qz_highlevel` | `geigen` | complete |
| `H` | `qz_highlevel` | `h_transpose` | complete |

Detailed source mappings and compatibility notes are recorded in `fpm.toml` and `docs/API_COVERAGE.md`.

## Compatibility notes

- R type/S3 dispatch is represented by Fortran generic interfaces selected at compile time.
- R list fields/classes and print methods are replaced by explicit derived result types.
- Real LAPACK packed eigenvectors are reconstructed into complex vectors using the same conjugate-pair convention as upstream QZ.
- Infinite/undefined generalized eigenvalues use IEEE infinity/NaN semantics so `ordqz` finite-only selectors behave like the R computation `ALPHA/BETA`.
- Eigenvector signs/phases and Schur-vector bases are not unique; tests therefore validate decomposition identities and eigenvector residuals rather than exact columns.

## Provenance and licensing

See `LICENSE`, `LICENSES/LAPACK_LICENSE.txt`, and `NOTICE.md`. The upstream R sources used for the translation are retained under `upstream/R/` for traceability. No dependency source is copied into this package.
