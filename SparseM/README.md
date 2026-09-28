# SparseM Fortran translation

This directory is a modern free-form Fortran translation of the computational core of the R package **SparseM 1.84-2** by Roger Koenker, with contributions credited upstream to Pin Tian Ng, Yousef Saad, Ben Shaby, and Martin Maechler. The upstream package is licensed under GPL (>= 2); this translation is distributed under **GPL-2.0-or-later** and keeps an upstream source snapshot for provenance.

The public Fortran API uses explicit sparse derived types for CSR, CSC, COO, symmetric sparse-row (SSR), and symmetric sparse-column (SSC) storage. CSR transpose and matrix products are sparse-aware. The default Cholesky path is now a maintained free-form translation of the central Ng-Peyton machinery used by SparseM: multiple-minimum-degree ordering, elimination-tree postordering, column counts, fundamental-supernode discovery, compressed supernodal symbolic factorization, numerical input, cache-block splitting, a level-1 reference update plus paired-target level-2/4/8 numerical kernels, and block triangular solves. Level 8 is the default, matching current SparseM. The factor retains the upstream-style `LINDX`/`XLINDX`/`LNZ`/`XLNZ` metadata and also exposes equivalent sparse CSR lower/upper factors for the rest of this Fortran API. `ordering="natural"` remains available as the earlier row-sparse fallback. No n-by-n dense factor is allocated. Several elementwise operations, binding/subsetting paths, and summary covariance calculations still use dense work arrays where that keeps the translation clear and deterministic.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated-package source is vendored. This translation is dependency-free, so there is no system `-lblas`/`-llapack` requirement and no copied shared dependency.

## Build and test

```text
fpm build
fpm test
```

The package is intended to build with gfortran on Windows as a top-level FPM project directory. `test/test_sparsem.f90` exercises the public computational API, while `test/test_ng_reference.f90` locks a compact Ng-Peyton factorization against upstream-reference ordering, symbolic storage, cache split, workspace size, and numerical `LNZ` values. `example/basic_sparsem.f90` demonstrates conversion and sparse linear least squares.

## Translation coverage

Package status: **substantial**. Coverage is **45 of 45 (100.0%)** exported computational R functions. The fraction measures computational functions with a meaningful Fortran mapping; it does **not** mean complete R/S4 compatibility or performance parity.

All exported computational R functions in the stated coverage basis now have a meaningful Fortran mapping. Plotting, printing, display methods, and pure class predicates are excluded from the coverage denominator under the stated coverage convention. A 100% mapping fraction does **not** imply exact R interface compatibility or sparse-factorization performance parity.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `as.matrix.csr` | `sparsem_core`, `sparsem_factor` | `csr_from_dense`, `csr_from_csc`, `csr_from_coo`, `csr_from_ssr`, `csr_from_ssc`, `chol_as_csr` | substantial |
| `as.matrix.csc` | `sparsem_core` | `csc_from_dense`, `csc_from_csr` | substantial |
| `as.matrix.ssr` | `sparsem_core` | `ssr_from_dense`, `ssr_from_csr` | substantial |
| `as.matrix.ssc` | `sparsem_core` | `ssc_from_dense`, `ssc_from_csr` | substantial |
| `as.matrix.coo` | `sparsem_core` | `coo_from_dense`, `coo_from_csr` | substantial |
| `as.matrix` | `sparsem_core` | `csr_to_dense`, `csc_to_dense`, `coo_to_dense`, `ssr_to_dense`, `ssc_to_dense` | substantial |
| `det` | `sparsem_factor` | `csr_determinant`, `chol_determinant` | substantial |
| `forwardsolve` | `sparsem_factor` | `chol_forwardsolve`, `chol_forwardsolve_matrix`, `chol_forwardsolve_csr` | substantial |
| `norm` | `sparsem_core` | `csr_norm` | complete |
| `rbind.matrix.csr` | `sparsem_core` | `csr_rbind` | substantial |
| `cbind.matrix.csr` | `sparsem_core` | `csr_cbind` | substantial |
| `Ops.matrix.csr` | `sparsem_core` | `csr_add`, `csr_subtract`, `csr_hadamard`, `csr_divide`, `csr_power`, `csr_intdiv`, `csr_modulo`, `csr_relation`, `csr_complement` | substantial |
| `Ops.matrix.diag.csr` | `sparsem_core` | `csr_add`, `csr_subtract`, `csr_hadamard`, `csr_divide`, `csr_power`, `csr_relation` | partial |
| `[.matrix.csr` | `sparsem_core` | `csr_submatrix` | substantial |
| `[<-.matrix.csr` | `sparsem_core` | `csr_assign_submatrix` | substantial |
| `[.matrix.diag.csr` | `sparsem_core` | `csr_submatrix` | partial |
| `slm` | `sparsem_factor` | `slm` | partial |
| `slm.fit` | `sparsem_factor` | `slm_fit` | substantial |
| `slm.wfit` | `sparsem_factor` | `slm_wfit` | substantial |
| `slm.fit.csr` | `sparsem_factor` | `slm_fit_csr`, `mslm_fit_csr` | substantial |
| `t` | `sparsem_core` | `csr_transpose` | substantial |
| `chol` | `sparsem_factor` | `csr_chol` | substantial |
| `coef.slm` | `sparsem_factor` | `slm_coef` | complete |
| `fitted.slm` | `sparsem_factor` | `slm_fitted` | complete |
| `residuals.slm` | `sparsem_factor` | `slm_residuals` | complete |
| `deviance.slm` | `sparsem_factor` | `slm_deviance` | substantial |
| `extractAIC.slm` | `sparsem_factor` | `slm_extract_aic` | substantial |
| `summary.slm` | `sparsem_factor` | `slm_summary` | substantial |
| `[<-.matrix.coo` | `sparsem_core` | `coo_assign_submatrix` | substantial |
| `[.matrix.coo` | `sparsem_core` | `coo_submatrix` | substantial |
| `dim` | `sparsem_core` | `csr_shape` | complete |
| `diff` | `sparsem_core` | `csr_diff_rows` | substantial |
| `diag` | `sparsem_core` | `csr_diag` | substantial |
| `diag<-` | `sparsem_core` | `csr_set_diag` | substantial |
| `determinant` | `sparsem_factor` | `csr_determinant`, `chol_determinant` | substantial |
| `backsolve` | `sparsem_factor` | `chol_backsolve`, `chol_backsolve_matrix`, `chol_backsolve_csr`, `chol_solve`, `chol_solve_matrix`, `chol_solve_csr` | substantial |
| `solve` | `sparsem_factor` | `csr_solve`, `csr_solve_matrix`, `csr_solve_csr_rhs`, `csr_inverse` | substantial |
| `%*%` | `sparsem_core` | `csr_matmul`, `csr_matvec`, `csr_vecmat` | substantial |
| `kronecker` | `sparsem_core` | `csr_kronecker` | substantial |
| `read.matrix.hb` | `sparsem_hb` | `read_matrix_hb` | substantial |
| `model.matrix` | `sparsem_hb` | `hb_model_matrix` | substantial |
| `model.response` | `sparsem_hb` | `hb_model_response` | substantial |
| `model.guess` | `sparsem_hb` | `hb_model_guess` | substantial |
| `model.xexact` | `sparsem_hb` | `hb_model_xexact` | substantial |
| `summary.mslm` | `sparsem_factor` | `mslm_summary` | substantial |

## Important compatibility notes

- `slm` accepts an explicit CSR model matrix and response vector. R formula parsing, model frames, contrasts, S3/S4 metadata, names, and NA actions are not translated.
- `chol`, `solve`, `det`, and `determinant` target symmetric positive-definite matrices, as does SparseM's `matrix.csr` solve path. The default `csr_chol` path translates the Ng-Peyton/Liu ordering, SFINIT-style elimination-tree analysis, SYMFCT-style compressed supernodal structure, INPNV/BFINIT-style numerical setup and cache splitting, block factorization, and BLKSLF/BLKSLB/BLKSLV-style triangular solves. `level=8` is the default; levels 1, 2, 4, and 8, `cache_kb`, `tiny`, `large`, `eps`, and optional `nnzlmax`/`nsubmax`/`tmpmax` ceilings are exposed. When ceilings are omitted, the R defaults are reproduced: `nsubmax=nnz(A)`, `nnzlmax=max(4*nnz(A), floor(0.2*nnz(A)^1.3))`, and `tmpmax=50*n`. `tmpmax` is checked before allocating the real update workspace. Tiny-pivot replacement returns status `16+k` while leaving a usable factor and recording `tiny_pivots=k`. Development oracle tests against the supplied fixed-form routines now include end-to-end driver comparisons on tridiagonal, star, grid, and irregular sparse SPD matrices: permutation, symbolic structures, cache splits, and final `LNZ` matched exactly at levels 1, 4, and 8. Separate `FNSPLT`/`FNTSIZ` checks matched cache partitioning and workspace size, including a 200-column supernode forced into multiple cache blocks. The level-2/4/8 matrix-update kernels follow the upstream paired-target scheduling, with explicit source-column unrolling; the modern level-2 kernel is numerically equivalent to the other modern levels but does not intentionally reproduce divergent behavior observed in the supplied legacy level-2 path. `chol_statistics` provides the factor nonzero/subscript counts, largest supernode sizes, temporary-workspace requirement, and the `LSTATS`-style factorization/triangular-solve operation counts. Cholesky objects cache both matrix determinant metadata and SparseM-compatible triangular-factor determinant metadata; `chol_as_csr` reproduces the upstream factor-to-CSR layout, and no-RHS `solve` is represented by `csr_inverse`. Dense and CSR right-hand sides are supported. R warning emission, elapsed-time fields, and S4 object construction are not reproduced. `csr_chol(..., ordering="natural")` retains the earlier row-sparse fallback.
- Harwell-Boeing input now supports the same central real/assembled use case as the upstream R reader: rectangular/unsymmetric or symmetric sparse coefficient storage, fixed-width integer/real fields, full RHS storage, and optional guess/exact-solution blocks. R connection objects and unsupported Harwell-Boeing matrix families are not reproduced.
- Multi-response fitting reuses one normal-equation Cholesky factor across response columns. `slm_summary` and `mslm_summary` now include coefficient standard errors, t statistics, two-sided p-values, coefficient correlations, R-squared measures, and F tests. R response naming, formula/effects metadata, weighted-model metadata, and list/class behavior remain omitted.
- Sparse division, exponentiation, integer division, modulo, some subsetting operations, and some binding paths use dense temporary arrays to preserve simple numerical semantics. They should not be treated as large-scale sparse implementations yet.
- The Fortran zero matrix uses zero stored entries. The R implementation's historical dummy-entry conventions are not reproduced.
- RNG behavior is irrelevant to the translated computational path; no random numbers are used to encode NaN/Inf sentinels.

## Provenance

See `NOTICE.md` and `upstream/`. The original package metadata and R source are retained solely for attribution, coverage accounting, and future parity work. The maintained Fortran source is entirely free-form `.f90`; no upstream fixed-form Cholesky source is copied into the package.
