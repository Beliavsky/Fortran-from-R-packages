# geigen - modern Fortran translation

This directory translates the computational API of the R package **geigen 2.4** to modern free-form Fortran. The upstream package computes generalized eigenvalues/eigenvectors, generalized Schur (QZ) decompositions, and generalized singular value decompositions (GSVDs) of matrix pairs.

## Build

The package is intended to live as a top-level directory in `Fortran-from-R-packages` and uses the FPM-compatible Beliavsky fork of Federico Perini's `fortran-lapack`:

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

The manifest does not link a separately installed BLAS or LAPACK library and no dependency source is copied into this package.

## Public Fortran API

Use `geigen_api`. It exports the common `dp` kind and the following typed numerical interfaces:

- `geigen(A,B,result[,symmetric,only_values])` for real or complex generalized eigensystems;
- `gqz(A,B,result[,sort])` for real or complex generalized Schur decompositions;
- `gevalues(result)` for generalized eigenvalues from a QZ result;
- `gsvd(A,B,result)` for real or complex rectangular matrix pairs;
- `gsvd_r`, `gsvd_or`, `gsvd_d1`, and `gsvd_d2` for the GSVD factors used by the upstream helper functions.

Results use derived types rather than R lists/S3 classes. Numerical or dimension failures are returned through each result's `info` field.

### GSVD implementation note

The current `fortran-lapack` API does not publicly export a `ggsvd3` generic. The translation therefore does not vendor the upstream fixed-form `D/ZGGSVD3` and `D/ZGGSVP3` sources. Instead it constructs the same GSVD identities from LAPACK operations that the dependency does expose: an economy SVD of the stacked matrix, a symmetric/Hermitian eigendecomposition of the cosine subspace, and an RQ factorization for the shared triangular factor. Tests validate the defining real and complex GSVD identities, including wide matrix pairs.

## Translation coverage

Package status: **substantial**.

**9 of 9 (100.0%)** exported computational R functions have meaningful Fortran mappings. The percentage measures mapped computational functions, not complete R-interface compatibility. The exported `gevalues` S3 dispatcher itself is not counted because it performs no numerical work; its two registered computational methods are counted. Printing methods, R S3 classes, `.Fortran` plumbing, and error-message formatting are outside the coverage basis.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `geigen` | `geigen_api` | `geigen` | substantial |
| `gqz` | `geigen_api` | `gqz` | substantial |
| `gevalues.xdgges` | `geigen_api` | `gevalues` | complete |
| `gevalues.xzgges` | `geigen_api` | `gevalues` | complete |
| `gsvd` | `geigen_api` | `gsvd` | substantial |
| `gsvd.R` | `geigen_api` | `gsvd_r` | complete |
| `gsvd.oR` | `geigen_api` | `gsvd_or` | complete |
| `gsvd.D1` | `geigen_api` | `gsvd_d1` | complete |
| `gsvd.D2` | `geigen_api` | `gsvd_d2` | complete |

Detailed source mappings and compatibility notes are recorded in `fpm.toml` and `docs/API_COVERAGE.md`.

## Compatibility notes

- R numeric/complex runtime dispatch is represented by Fortran generic resolution.
- Symmetric/Hermitian detection is numerical and follows the upstream automatic choice unless the caller overrides it.
- Real LAPACK packed eigenvectors are reconstructed into complex conjugate-pair columns.
- Generalized eigenvalues retain the upstream special handling for zero and non-finite beta values.
- `gqz` supports the upstream ordering codes `N`, `-`, `+`, `S`, `B`, and `R`.
- Eigenvector signs/phases, Schur bases, and GSVD bases are not unique; tests validate residual/decomposition identities rather than exact columns.
- The translated GSVD uses supported LAPACK primitives rather than `GGSVD3`; small ordering/sign differences relative to the upstream LAPACK routine are therefore possible while the defining decompositions remain equivalent.

## Provenance and licensing

See `LICENSE`, `NOTICE.md`, `NOTICE_UPSTREAM.txt`, and `LICENSES/LAPACK_LICENSE.txt`. Upstream R sources are retained under `upstream/R/` for traceability. The fixed-form LAPACK sources distributed by the R package are intentionally not copied into this translation.
