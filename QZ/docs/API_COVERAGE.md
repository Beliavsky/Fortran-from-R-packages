# API coverage

Coverage basis: exported computational R functions.

Package status: **substantial**.

Mapped: **18 of 18 (100.0%)**.

The coverage fraction records whether an exported computational operation has a meaningful Fortran implementation. It does not claim identical R objects, S3 classes, error messages, dynamic-library plumbing, workspace-query conventions, or print behavior.

| R function | Fortran mapping | Status | Main compatibility note |
|---|---|---|---|
| `qz.dgees` | `qz_lapack:qz_dgees` | complete | R list/S3 wrapper omitted. |
| `qz.dgeev` | `qz_lapack:qz_dgeev` | complete | Complex vectors reconstructed from LAPACK real pair storage. |
| `qz.dgges` | `qz_lapack:qz_dgges` | complete | R list/S3 wrapper omitted. |
| `qz.dggev` | `qz_lapack:qz_dggev` | complete | Complex vectors reconstructed from LAPACK real pair storage. |
| `qz.dtgsen` | `qz_lapack:qz_dtgsen` | complete | Uses explicit TGSEN interface. |
| `qz.dtrsen` | `qz_lapack:qz_dtrsen` | complete | Uses explicit TRSEN interface. |
| `qz.zgees` | `qz_lapack:qz_zgees` | complete | R list/S3 wrapper omitted. |
| `qz.zgeev` | `qz_lapack:qz_zgeev` | complete | Uses explicit complex GEEV interface. |
| `qz.zgges` | `qz_lapack:qz_zgges` | complete | Uses explicit complex GGES interface. |
| `qz.zggev` | `qz_lapack:qz_zggev` | complete | Uses explicit complex GGEV interface. |
| `qz.ztgsen` | `qz_lapack:qz_ztgsen` | complete | Uses explicit complex TGSEN interface. |
| `qz.ztrsen` | `qz_lapack:qz_ztrsen` | complete | Uses explicit complex TRSEN interface. |
| `fda.geigen` | `qz_highlevel:fda_geigen` | substantial | Numerical augmented problem and normalization retained; R objects/errors differ. |
| `ordqz` | `qz_highlevel:ordqz` | substantial | All upstream keywords plus explicit cluster masks; Fortran generic dispatch. |
| `qz` | `qz_highlevel:qz_decompose` | substantial | Standard/generalized and real/complex dispatch via generics. |
| `qz.geigen` | `qz_highlevel:qz_geigen` | substantial | Standard/generalized and real/complex dispatch via generics. |
| `geigen` | `qz_highlevel:geigen` | complete | Alias generic for `qz_geigen`. |
| `H` | `qz_highlevel:h_transpose` | complete | Direct Hermitian transpose. |

Excluded from the denominator: `print.*` S3 methods and `.onLoad`/`.onUnload` runtime hooks.
