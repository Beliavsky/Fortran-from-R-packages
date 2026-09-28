# API coverage

## Counting basis

Coverage follows the project convention: distinct exported computational R functions, including registered computational S3 methods and excluding constants, presentation-only methods, destructor/memory-management helpers, and the internal `.get_cnst` helper.

`NAMESPACE` registers six computational S3 methods in addition to the constructor and exported `Rmumps__*` wrappers: `solve.Rcpp_Rmumps`, `solvet.Rcpp_Rmumps`, `determinant.Rcpp_Rmumps`, `dim.Rcpp_Rmumps`, `nrow.Rcpp_Rmumps`, and `ncol.Rcpp_Rmumps`. These are included in the denominator as required by the translation-coverage convention.

Package status: **partial**.

Mapped exported computational functions: **13 of 13 (100%)**. The percentage counts mapped functions, not complete R or MUMPS compatibility.

| R function | Fortran implementation | Status | Main difference |
|---|---|---|---|
| `Rmumps` | `rmumps_init_triplet`, `rmumps_init_triplet_zero_based`, and `rmumps_matrix` | partial | Triplet construction and copy-state option; no Rcpp/S4 conversion or MUMPS sparse backend |
| `Rmumps__solveptr` | `rmumps_solve_vector_inplace`, `rmumps_solve_matrix_inplace` | substantial | Fortran arrays replace external pointers |
| `Rmumps__ptr_ijv` | `rmumps_init_triplet_zero_based`, `rmumps_init_triplet` | substantial | Zero-based upstream indexing and one-based Fortran indexing are both supported; no pointer ownership layer |
| `Rmumps__triplet` | `rmumps_get_triplet` | substantial | Allocatable arrays replace an R `simple_triplet_matrix` object |
| `Rmumps__set_mat_ptr` | `rmumps_set_values` | substantial | Copies values instead of aliasing a raw pointer |
| `Rmumps__set_permutation` | `rmumps_set_permutation` | partial | Ordering code stored/validated but dense LU does not apply MUMPS ordering algorithms |
| `Rmumps__get_permutation` | `rmumps_get_permutation` | complete | Direct state query |
| `solve.Rcpp_Rmumps` | dense/sparse `rmumps_solve_*`, `rmumps_inverse` | substantial | Dense backend; no R Matrix/S4 coercion or copy flag semantics |
| `solvet.Rcpp_Rmumps` | dense/sparse `rmumps_solvet_*`, `rmumps_inverse_transpose` | substantial | Dense backend; same cached LU is reused for transpose systems |
| `determinant.Rcpp_Rmumps` | `rmumps_determinant` | substantial | Computes sign and modulus directly; no R determinant-class object |
| `dim.Rcpp_Rmumps` | `rmumps_dim` | complete | Direct shape query |
| `nrow.Rcpp_Rmumps` | `rmumps_nrow` | complete | Direct row-count query |
| `ncol.Rcpp_Rmumps` | `rmumps_ncol` | complete | Direct column-count query |

## Additional translated class behavior

The module also exposes computational behavior corresponding to important public Rcpp class methods that are not separate namespace-level R functions in the counting denominator:

- `symbolic()` -> `rmumps_symbolic`
- `numeric()` -> `rmumps_numeric`
- `solvev()` / `solvem()` -> dense `rmumps_solve_*`
- `solves()` / `solvestm()` -> `rmumps_solve_sparse_rhs`
- `inv()` -> `rmumps_inverse`
- `det()` -> `rmumps_det`
- `log2det()` -> `rmumps_log2det` with an explicit determinant-sign output
- `set_mat_data()` -> `rmumps_set_values`
- `triplet()` -> `rmumps_get_triplet`
- symmetry field query -> `rmumps_get_sym`
- `set_icntl()` / `get_icntl()` -> `rmumps_set_icntl` / `rmumps_get_icntl`
- `set_cntl()` / `get_cntl()` -> `rmumps_set_cntl` / `rmumps_get_cntl`
- `set_keep()` / `get_keep()` -> `rmumps_set_keep` / `rmumps_get_keep`
- `get_infos()` -> `rmumps_get_infos`
- `rhs` / `mrhs` properties -> `rmumps_set_rhs`, `rmumps_get_rhs`, `rmumps_set_mrhs`, `rmumps_get_mrhs`
- `copy` field -> `rmumps_set_copy`, `rmumps_get_copy`
- `ncore` field -> `rmumps_set_ncore`, `rmumps_get_ncore`
- `set_perm_in()` -> `rmumps_set_perm_in`

The control arrays are initialized with the explicit defaults established by upstream `tri_init`. The dense fallback stores changes to these arrays, but only explicit high-level operations with a meaningful fallback analogue alter the numerical path. Direct `set_icntl()` calls preserve cached factors just as the upstream wrapper does; the higher-level `set_permutation()` method invalidates them when its ordering code changes. `ICNTL(9)=2` is honored by ordinary dense/sparse solve paths as a transpose solve request, and RHS-mode state in `ICNTL(20)` follows dense versus sparse calls. `rmumps_get_infos` returns the same array shapes as MUMPS and populates the available status plus determinant mantissa/exponent fields; unsupported MUMPS diagnostics remain zero. The `copy` flag cannot reproduce R aliasing in an allocatable-array Fortran API, and `ncore` is stored but does not parallelize the serial dense fallback. Successful ordinary solves update the stored `rhs` or `mrhs` state. `rmumps_set_perm_in` validates and stores a one-based user permutation, sets `ICNTL(7)=1` and `ICNTL(28)=1`, and invalidates analysis/factors, but the dense LU backend does not apply that permutation.

Sparse RHS support accepts triplets rather than R `dgCMatrix` or `simple_triplet_matrix` objects. Duplicate RHS triplets are summed. A `0 x 0` sparse RHS requests an inverse, reproducing the upstream special case. Singular dense-fallback factorizations are surfaced as MUMPS-compatible status `-10`. The KEEP defaults exercised by the upstream tests (`KEEP(10)` and `KEEP(34)`) are initialized from the compiler default-integer storage size.

## Major remaining compatibility gaps

The following public Rcpp-module features remain MUMPS-specific or R-interface-specific and are not implemented as compatible Fortran operations yet:

- MUMPS sparse multifrontal symbolic/numeric factorization and memory scaling;
- AMD/AMF/SCOTCH/PORD/METIS/QAMD ordering algorithms;
- numerical effects of most MUMPS `ICNTL`, `CNTL`, and `KEEP` controls, and MUMPS-only diagnostic fields in `INFO`, `RINFO`, `INFOG`, and `RINFOG`;
- applying user-supplied MUMPS permutation vectors during factorization and returning MUMPS-generated symbolic/unsymmetric permutations;
- actual OpenMP use of the stored core-count control and distributed-memory behavior;
- exact `copy=FALSE` R object/pointer aliasing semantics;
- Rcpp/S4/Matrix/slam construction, coercion, and object-return behavior;
- reporting an actual MUMPS library version, because the current backend does not link MUMPS.

These differences are why package status remains `partial` despite 100% mapped-function coverage.
