# matrixcalc

Modern free-form Fortran translation of the computational routines in the R package **matrixcalc** 1.0-6 by Frederick Novomestky. The upstream package is licensed GPL (>= 2); this translation is distributed under GPL-2.0-or-later and retains the upstream R sources under `upstream/` for provenance.

The implementation provides matrix transformations (`vec`, `vech`, duplication/elimination/commutation matrices), direct/Kronecker products and sums, matrix norms and predicates, LU/inverse/power/rank operations, and the package's special-matrix constructors. Linear solves and symmetric eigendecompositions reuse the sibling `rfortran-linalg` package; no BLAS or LAPACK source is vendored here and the manifest does not link system BLAS/LAPACK.

## Build

Place this directory at the repository root alongside `rfortran-core` and `rfortran-linalg`, then run:

```text
fpm build
fpm test
fpm run --example basic
```

The public facade is `use matrixcalc`. Real calculations consistently use `real(dp)`, with `dp` re-exported from `r_kinds`.

## Translation coverage

**Package status:** substantial  
**Coverage:** **60 of 60 (100.0%)** exported computational R functions.

The percentage measures mapped computational functions, not complete R compatibility. In particular, Fortran uses typed arrays/derived types instead of R lists, reports selected numerical failures through status values, and does not reproduce R's dynamic type validation or exact error strings.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `D.matrix` | `matrixcalc_structured` | `d_matrix` | complete |
| `E.matrices` | `matrixcalc_structured` | `e_matrices` | substantial |
| `H.matrices` | `matrixcalc_structured` | `h_matrices` | substantial |
| `K.matrix` | `matrixcalc_structured` | `k_matrix` | complete |
| `L.matrix` | `matrixcalc_structured` | `l_matrix` | complete |
| `N.matrix` | `matrixcalc_structured` | `n_matrix` | complete |
| `T.matrices` | `matrixcalc_structured` | `t_matrices` | substantial |
| `commutation.matrix` | `matrixcalc_structured` | `commutation_matrix` | complete |
| `creation.matrix` | `matrixcalc_special` | `creation_matrix` | complete |
| `direct.prod` | `matrixcalc_core` | `direct_prod` | complete |
| `direct.sum` | `matrixcalc_core` | `direct_sum` | complete |
| `duplication.matrix` | `matrixcalc_structured` | `duplication_matrix` | substantial |
| `elimination.matrix` | `matrixcalc_structured` | `elimination_matrix` | complete |
| `entrywise.norm` | `matrixcalc_core` | `entrywise_norm` | substantial |
| `fibonacci.matrix` | `matrixcalc_special` | `fibonacci_matrix` | complete |
| `frobenius.matrix` | `matrixcalc_special` | `frobenius_matrix` | complete |
| `frobenius.norm` | `matrixcalc_core` | `frobenius_norm` | complete |
| `frobenius.prod` | `matrixcalc_core` | `frobenius_prod` | complete |
| `hadamard.prod` | `matrixcalc_core` | `hadamard_prod` | complete |
| `hankel.matrix` | `matrixcalc_special` | `hankel_matrix` | complete |
| `hilbert.matrix` | `matrixcalc_special` | `hilbert_matrix` | complete |
| `hilbert.schmidt.norm` | `matrixcalc_core` | `hilbert_schmidt_norm` | complete |
| `inf.norm` | `matrixcalc_core` | `inf_norm` | complete |
| `is.diagonal.matrix` | `matrixcalc_linalg` | `is_diagonal_matrix` | complete |
| `is.idempotent.matrix` | `matrixcalc_linalg` | `is_idempotent_matrix` | complete |
| `is.indefinite` | `matrixcalc_linalg` | `is_indefinite` | substantial |
| `is.negative.definite` | `matrixcalc_linalg` | `is_negative_definite` | substantial |
| `is.negative.semi.definite` | `matrixcalc_linalg` | `is_negative_semi_definite` | substantial |
| `is.non.singular.matrix` | `matrixcalc_linalg` | `is_non_singular_matrix` | substantial |
| `is.positive.definite` | `matrixcalc_linalg` | `is_positive_definite` | substantial |
| `is.positive.semi.definite` | `matrixcalc_linalg` | `is_positive_semi_definite` | substantial |
| `is.singular.matrix` | `matrixcalc_linalg` | `is_singular_matrix` | substantial |
| `is.skew.symmetric.matrix` | `matrixcalc_linalg` | `is_skew_symmetric_matrix` | complete |
| `is.square.matrix` | `matrixcalc_linalg` | `is_square_matrix` | complete |
| `is.symmetric.matrix` | `matrixcalc_linalg` | `is_symmetric_matrix` | complete |
| `lower.triangle` | `matrixcalc_core` | `lower_triangle` | complete |
| `lu.decomposition` | `matrixcalc_linalg` | `lu_decomposition` | substantial |
| `matrix.inverse` | `matrixcalc_linalg` | `matrix_inverse` | substantial |
| `matrix.power` | `matrixcalc_linalg` | `matrix_power` | substantial |
| `matrix.rank` | `matrixcalc_linalg` | `matrix_rank` | substantial |
| `matrix.trace` | `matrixcalc_core` | `matrix_trace` | complete |
| `maximum.norm` | `matrixcalc_core` | `maximum_norm` | substantial |
| `one.norm` | `matrixcalc_core` | `one_norm` | complete |
| `pascal.matrix` | `matrixcalc_special` | `pascal_matrix` | substantial |
| `%s%` | `matrixcalc_core` | `direct_sum` | complete |
| `set.submatrix` | `matrixcalc_core` | `set_submatrix` | substantial |
| `shift.down` | `matrixcalc_core` | `shift_down` | complete |
| `shift.left` | `matrixcalc_core` | `shift_left` | complete |
| `shift.right` | `matrixcalc_core` | `shift_right` | complete |
| `shift.up` | `matrixcalc_core` | `shift_up` | complete |
| `spectral.norm` | `matrixcalc_linalg` | `spectral_norm` | substantial |
| `stirling.matrix` | `matrixcalc_special` | `stirling_matrix` | substantial |
| `svd.inverse` | `matrixcalc_linalg` | `svd_inverse` | substantial |
| `symmetric.pascal.matrix` | `matrixcalc_special` | `symmetric_pascal_matrix` | complete |
| `toeplitz.matrix` | `matrixcalc_special` | `toeplitz_matrix` | complete |
| `u.vectors` | `matrixcalc_structured` | `u_vectors` | substantial |
| `upper.triangle` | `matrixcalc_core` | `upper_triangle` | complete |
| `vandermonde.matrix` | `matrixcalc_special` | `vandermonde_matrix` | complete |
| `vec` | `matrixcalc_core` | `vec` | substantial |
| `vech` | `matrixcalc_core` | `vech` | substantial |

Detailed notes for every mapping are in `fpm.toml` and `docs/API_COVERAGE.md`.
