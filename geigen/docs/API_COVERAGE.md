# API coverage

Coverage basis: distinct exported computational R functions. Non-computational S3 printing and the `gevalues` dispatcher are excluded; the two registered computational `gevalues` methods are included.

Package status: **substantial**.

Mapped: **9 of 9 (100.0%)**.

| R function | R source | Fortran procedure | Status | Compatibility note |
|---|---|---|---|---|
| `geigen` | `upstream/R/geigen.R` | `geigen_api::geigen` | substantial | Real/complex and symmetric/Hermitian paths translated; R coercion/list behavior omitted. |
| `gqz` | `upstream/R/gqz.R` | `geigen_api::gqz` | substantial | QZ factors and sort codes translated; S3 result classes omitted. |
| `gevalues.xdgges` | `upstream/R/gqz.R` | `geigen_api::gevalues` | complete | Typed QZ result replaces R S3 dispatch. |
| `gevalues.xzgges` | `upstream/R/gqz.R` | `geigen_api::gevalues` | complete | Typed QZ result replaces R S3 dispatch. |
| `gsvd` | `upstream/R/gsvd.R` | `geigen_api::gsvd` | substantial | GSVD constructed from supported LAPACK primitives because `la_lapack` does not export `ggsvd3`; defining identities are retained. |
| `gsvd.R` | `upstream/R/gsvd.R` | `geigen_api::gsvd_r` | complete | Returns common triangular factor. |
| `gsvd.oR` | `upstream/R/gsvd.R` | `geigen_api::gsvd_or` | complete | Returns `[0 R]`. |
| `gsvd.D1` | `upstream/R/gsvd.R` | `geigen_api::gsvd_d1` | complete | Returns first GSVD diagonal/block factor. |
| `gsvd.D2` | `upstream/R/gsvd.R` | `geigen_api::gsvd_d2` | complete | Returns second GSVD diagonal/block factor. |

The coverage fraction is a function-mapping metric and does not mean complete compatibility with R object classes, coercion, diagnostics, or printing.
