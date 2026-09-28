# rARPACK

Modern free-form Fortran implementation of the numerical API documented by the
R package `rARPACK` 0.11-0, packaged with FPM.

The upstream R release is a compatibility shell around `RSpectra`: its three R
exports contain no eigensolver implementation of their own. This translation
therefore keeps the upstream wrapper sources for provenance while implementing
useful Fortran counterparts with the repository's shared `rfortran-arpack`
backend. ARPACK, BLAS, LAPACK, `r.f90`, `r_mod.f90`, and translated dependency
sources are not vendored here.

## Build

Place this directory at repository top level beside `rfortran-arpack`, then run:

```text
fpm build
fpm test
fpm run --example basic
```

The manifest uses:

```toml
[dependencies]
rfortran-arpack = { path = "../rfortran-arpack" }
```

No system `-larpack`, `-llapack`, or `-lblas` link flags are used by this
package.

## Public API

```fortran
use rarpack, only: dp, eigs, eigs_result, eigs_sym, eigs_sym_result, svds, svds_result
```

`eigs` computes selected eigenpairs of a dense real nonsymmetric matrix and
returns complex values/vectors. `eigs_sym` handles dense real symmetric
problems and can read either triangle. Both provide real shift-and-invert mode.
`svds` computes leading singular triplets of a dense real rectangular matrix.

Optional `ncv`, `tol`, and `maxitr` arguments correspond to the main numerical
controls in the R `opts` list. `eigs` and `eigs_sym` additionally accept
`retvec`; `eigs_sym` accepts `lower`; `svds` accepts `nu` and `nv`.

The maintained Fortran code uses one kind, `dp = real64`, and deterministic
starting residuals. No fast-math assumptions are required or enabled.

## Translation coverage

Package status: **substantial**.

Strict translation coverage: **0 of 0 (0%)** exported computational R
functions. The fraction measures mapped computational functions rather than
complete R compatibility. Under the requested convention, functions that
merely dispatch to another package are excluded; all three rARPACK 0.11-0
exports are direct RSpectra dispatchers. For this empty denominator, the
manifest records `0.0` as the numeric fraction convention.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| None under the strict coverage basis | - | - | - |

For practical use, the three documented operations have substantial Fortran
compatibility implementations even though they are excluded from the strict R
wrapper-function fraction:

| Upstream export | Fortran module | Fortran procedure(s) | Compatibility status |
| --- | --- | --- | --- |
| `eigs` | `rarpack` | `eigs` | substantial |
| `eigs_sym` | `rarpack` | `eigs_sym` | substantial |
| `svds` | `rarpack` | `svds` | substantial |

See `API_COVERAGE.md` for supported selectors, result types, and material
interface differences.

## Important compatibility differences

The Fortran API currently supports dense real matrices. R `Matrix` sparse
classes, R closures as implicit matrix operators, S3/R object construction,
warning behavior, and R list/recycling semantics are not reproduced. General
nonsymmetric eigenvalues and eigenvectors are returned as `complex(dp)`.
Shift-and-invert accepts a real scalar shift and factors the dense shifted
matrix internally with partial pivoting.

`svds` applies ARPACK to the smaller normal-equation operator without explicitly
forming that Gram matrix. This is appropriate for a compact compatibility
implementation but is not an exact reproduction of RSpectra's algorithms or
conditioning behavior.

## License and provenance

See `LICENSE`, `NOTICE.md`, and `upstream/`. Historical ARPACK and arpack-ng
copyright notices from the R package are preserved verbatim in
`upstream/inst/COPYRIGHTS`.
