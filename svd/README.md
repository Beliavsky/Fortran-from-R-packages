# svd

Modern free-form Fortran/FPM translation of the computational API of the R package
`svd` 0.5.8 by Anton Korobeynikov and contributors.

The upstream package exposes matrix-free multiplication plus leading singular-value
and eigenvalue solvers backed by PROPACK and nuTRLan. This translation intentionally
does **not** copy those bundled numerical libraries. It reuses the existing sibling
`RSpectra` translation for real dense/matrix-free iterative SVD/eigensolver work and
`rfortran-linalg` for the complex dense SVD. Both shared packages use the repository's
pinned pure-Fortran numerical backend, so this package has no system BLAS/LAPACK/ARPACK
link requirement in its FPM manifest.

## Build

Place this directory at the repository root beside `RSpectra`, `rfortran-linalg`, and
their existing sibling dependencies, then run:

```text
fpm build
fpm test
fpm run svd_demo
fpm run --example dense_svd
fpm run --example matrix_free_svd
```

The project uses only free-form maintained Fortran. No translated dependency source is
vendored here.

## Main API

- `extmat_callbacks`, `extmat_operator`, `make_extmat`: matrix-free real operators.
- `ematmul`, `materialize_extmat`, `extmat_nrow`, `extmat_ncol`, `is_extmat`.
- `propack_svd`: generic interface for dense or matrix-free real leading SVDs.
- `trlan_svd`: dense/operator singular values plus left singular vectors through the
  left normal operator `A*A^T`, matching the upstream result shape.
- `trlan_eigen`: dense/operator leading real symmetric eigenpairs.
- `ztrlan_svd_complex`: dense conventional complex SVD returning leading singular
  values and left singular vectors.

All public real arrays use the package-wide `dp = real64` kind defined once in
`svd_kinds`. This is the same real kind used by the sibling numerical backends.

## Algorithmic compatibility

The public numerical operations are preserved, but the iterative engine is deliberately
shared rather than copied. `propack.svd` therefore maps its Krylov dimension/tolerance/
restart controls onto the existing RSpectra/ARPACK interface instead of reproducing
PROPACK's partial-reorthogonalization implementation. Likewise, the real TRLan routines
solve the same operator problems with ARPACK rather than nuTRLan thick restart.

Upstream TRLan accepts previously computed `lambda`/`U` arrays as locked vectors. The
Fortran interfaces accept the same numerical information, validate its shape, and use
one matched vector to initialize the iterative backend; the full locked-subspace behavior
is not reproduced.

For complex `ztrlan.svd`, the translation uses a standard Hermitian complex SVD through
`rfortran-linalg`. The upstream native implementation applies an unconjugated transpose
inside its complex normal operator; reproducing that behavior would not implement the
usual mathematical complex SVD and is therefore intentionally not used here.

## Translation coverage

Package status: **substantial**.

**12 of 12 (100%)** exported computational R functions have meaningful Fortran mappings.
This percentage measures mapped computational functions, not complete compatibility with
R's S4 objects, environments, warnings, RNG, bundled solver internals, or all option-list
semantics. Presentation/coercion S4 methods are outside the coverage denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `is.extmat` | `svd_extmat` | `is_extmat` | complete |
| `extmat.ncol` | `svd_extmat` | `extmat_ncol` | complete |
| `extmat.nrow` | `svd_extmat` | `extmat_nrow` | complete |
| `ematmul` | `svd_extmat` | `ematmul` | complete |
| `extmat` | `svd_extmat` | `make_extmat` | substantial |
| `as.array.extmat` | `svd_extmat` | `materialize_extmat` | substantial |
| `as.matrix.extmat` | `svd_extmat` | `materialize_extmat` | substantial |
| `as.vector.extmat` | `svd_extmat` | `extmat_vector` | substantial |
| `propack.svd` | `svd_solvers` | `propack_svd_dense`, `propack_svd_operator` | substantial |
| `trlan.svd` | `svd_solvers` | `trlan_svd_dense`, `trlan_svd_operator` | substantial |
| `trlan.eigen` | `svd_solvers` | `trlan_eigen_dense`, `trlan_eigen_operator` | substantial |
| `ztrlan.svd` | `svd_solvers` | `ztrlan_svd_complex` | substantial |

The three explicitly registered S3 conversion methods are counted because they perform
materialization/vectorization. The non-exported upstream `ztrlan.eigen` function is not
counted. R-only S4 `exportMethods` dispatch (coercion, transpose, `%*%`, `crossprod`, and
`tcrossprod`) is excluded as interface machinery under the task's instruction to skip
R-specific interfaces; equivalent numerical work is available through the typed operator
API, `ematmul`, `materialize_extmat`, and ordinary Fortran `matmul`.

## License and provenance

See `NOTICE.md` and `provenance/`. The upstream DESCRIPTION declares BSD-3-Clause plus
its `LICENSE` metadata, while the R/C interface files being computationally reimplemented
carry GPL-2-or-later notices. The maintained Fortran compatibility layer conservatively
uses GPL-2.0-or-later. PROPACK and nuTRLan are **not** copied or linked into this package;
their original license texts are retained under `provenance/licenses/` solely to preserve
upstream attribution and provenance.
