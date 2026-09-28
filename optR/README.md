# optR - modern Fortran translation

This directory translates the computational core of the R package **optR 1.2.5** to modern free-form Fortran. The upstream package, by Prakash (PKS Prakash), solves linear systems with Gaussian elimination, LU decomposition, Gauss-Seidel iteration, conjugate-gradient-style iteration, and Cholesky decomposition, and also provides matrix inversion, projection matrices, finite-difference Jacobians, and Newton-Raphson root finding.

The maintained Fortran source is dependency-free. This is intentional: the numerical algorithms are the substance of optR itself, so replacing them with an external BLAS/LAPACK solver would remove the code being translated. No system BLAS/LAPACK link flags are used.

## Build and test

From the `optR` directory:

```text
fpm build
fpm test
fpm run --example solve_demo
```

The source uses a single real kind, `dp = real64`, defined in `optr_kinds` and re-exported by `optr`. It contains no fixed-form Fortran and requires no separately installed numerical library.

## Public Fortran API

Import the user-facing module with:

```fortran
use optr
```

The main entry points are:

- `optr_solve`: dispatches to `gauss`, `LU`, `gaussseidel`, `cgm`, or `choleski`/`cholesky` for vector or matrix right-hand sides.
- `optr_fit_design`: fits a numeric design matrix directly, including optional observation weights; rectangular problems follow upstream optR by solving normal equations.
- `gaussian_solve`, `lu_solve`, `cholesky_solve`: direct linear solvers.
- `lu_decompose`, `lu_split`, `cholesky_decompose`: matrix factorizations.
- `gauss_seidel`, `conjugate_gradient`: iterative solvers.
- `invert_matrix`, `hat_matrix`: matrix inverse and ordinary/generalized projection matrices.
- `finite_difference_jacobian`, `newton_raphson`: nonlinear-system tools using a `vector_function` callback.
- `predict_linear`: matrix-based prediction from fitted coefficients.

Procedures return integer status codes rather than raising R errors. `optr_success` is zero; other public status constants distinguish shape errors, singularity, non-positive-definiteness, nonconvergence, and invalid arguments.

## Numerical compatibility notes

The translation follows the upstream algorithms and defaults where they materially define the numerical operation. It deliberately makes several Fortran-oriented changes:

- Direct Gaussian and LU paths use explicit pivot checks and stable row pivoting rather than depending on R warnings or non-finite arithmetic.
- `inv.optR` is represented by `invert_matrix`, which factors the matrix once and solves all identity columns.
- The upstream Gauss-Seidel relaxation formula is retained; its `w` is a multiplicative factor and is not silently replaced by a different SOR definition.
- The upstream CGM recurrence is retained rather than being replaced by an unrelated library solver.
- R formulas, data frames, model frames, contrasts, S3 classes, plotting, printing, and interactive/presentation behavior are outside the Fortran API.
- Callback-based R functions used by `jacobian`/Newton-Raphson are represented by the `vector_function` procedure-pointer wrapper. The callback interface is pure so the numerical routines can themselves remain pure.

See `API_COVERAGE.md` for internal-kernel coverage and compatibility details. The complete uploaded R source snapshot is retained under `upstream/` for provenance and comparison.

## Translation coverage

Package status: **substantial**.

Coverage is **11 of 11 (100%)** exported computational R functions. This fraction measures functions with meaningful Fortran mappings, not complete R-interface compatibility. In particular, formula parsing and S3 object behavior are intentionally omitted even when the underlying numerical operation is mapped.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `LUsplit` | `optr` | `lu_split` | complete |
| `choleskiDecomposition` | `optr` | `cholesky_decompose` | complete |
| `hatMatrix` | `optr` | `hat_matrix` | substantial |
| `inv.optR` | `optr` | `invert_matrix` | substantial |
| `jacobian` | `optr` | `finite_difference_jacobian` | substantial |
| `optR` | `optr` | `optr_solve`, `optr_fit_design` | substantial |
| `optRFun` | `optr` | `newton_raphson` | substantial |
| `optRFun.newtonRapson` | `optr` | `newton_raphson` | substantial |
| `optR.default` | `optr` | `optr_solve`, `optr_fit_design` | substantial |
| `optR.formula` | `optr` | `optr_fit_design` | partial |
| `predict.optR` | `optr` | `predict_linear` | substantial |

`print.optR` and `summary.optR` are excluded from the denominator because they are presentation/plotting methods. Internal helpers are also excluded from the denominator, although their numerical work is translated where needed.

## License and provenance

The upstream package declares `GPL (>= 2)`. The Fortran translation is therefore distributed under **GPL-2.0-or-later**. See `LICENSE` and `NOTICE.md`. Upstream source and documentation are retained under `upstream/` without changing their authorship.
