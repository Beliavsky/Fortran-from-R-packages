# Translation notes

## Upstream

This is a standalone modern Fortran translation of the computational portions
of **psqn 0.3.2**, by Benjamin Christoffersen. The upstream package declares
`Apache License (>= 2)`. The derivative Fortran sources use the SPDX identifier
`Apache-2.0`; the complete license text is retained in `LICENSE` and the
upstream package metadata is retained in `original/DESCRIPTION`.

## Source mapping

| Upstream source | Fortran translation | Notes |
| --- | --- | --- |
| `inst/include/psqn-misc.h`, `constant.h` | `src/psqn_types.f90` | Status constants, options, callback interfaces, result types |
| `inst/include/lp.h` | `src/psqn_linalg.f90` | Packed symmetric algebra, Kahan summation, rank-one/BFGS updates, Cholesky utilities |
| `inst/include/intrapolate.h` | `src/psqn_interpolation.f90` | Line-search interpolation |
| `inst/include/richardson-extrapolation.h` | `src/psqn_richardson.f90` | Richardson-extrapolated vector derivatives |
| `inst/include/psqn-bfgs.h` | `src/psqn_bfgs.f90`, `src/psqn_linesearch.f90` | Dense inverse-BFGS and Wolfe line search |
| `inst/include/psqn.h` | `src/psqn_core.f90` | Partially separable quasi-Newton, BFGS/SR1 element updates, CG, masks, preconditioning, augmented Lagrangian, Hessians, private optimization |
| `src/r-api.cpp`, `src/RcppExports.cpp`, `R/RcppExports.R` | Not translated | R/Rcpp interface and registration rather than numerical algorithms |
| `inst/include/psqn-Rcpp-wrapper.h` | Not translated | Rcpp object adapter |
| `inst/include/psqn-reporter.h` | Not translated | Platform/reporting support |
| Vignettes and presentation material | Not translated | Non-computational; plotting/presentation intentionally skipped |

## API differences

* Fortran indices in `psqn_element_spec%idx` and `masked` are **1-based**.
* The port is standalone and has no R, Rcpp, RcppEigen, Eigen, or Matrix-package
  dependency.
* The upstream Eigen incomplete-Cholesky preconditioner is replaced by a
  standalone **regularized dense Cholesky** preconditioner. This preserves the
  preconditioned-CG role while eliminating the Eigen dependency.
* `psqn_pre_block` is implemented for the structured global/private API. In the
  generic API, where no global/private partition is known, a block request is
  handled as diagonal preconditioning.
* The current Fortran port is serial. Upstream optional OpenMP execution is an
  execution optimization, not a distinct numerical method, and is not included
  in this first standalone version.
* A very small directional-derivative guard is used by the Wolfe line search
  when roundoff leaves an already-stationary iterate with a nonzero derivative.
  This prevents spurious line-search failure at machine precision.

## Public computational entry points

The umbrella module `psqn` exports:

* `psqn_bfgs_optimize`
* `psqn_optimize_generic`
* `psqn_optimize_structured`
* `psqn_optimize_private_structured`
* `psqn_aug_lagrang_generic`
* `psqn_aug_lagrang_structured`
* `psqn_generic_hess`
* `psqn_structured_hess`
* `psqn_make_structured_specs`
* `richardson_vector_derivative`

The element callback interface is

```fortran
subroutine element_eval(i, x, f, g, comp_grad)
  integer, intent(in) :: i
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: f
  real(dp), intent(out) :: g(:)
  logical, intent(in) :: comp_grad
end subroutine element_eval
```

`i` identifies the element function, and `x`/`g` have the local dimension given
by the corresponding `psqn_element_spec%idx` entry.

## R API correspondence

| R routine | Fortran routine |
| --- | --- |
| `psqn()` | `psqn_optimize_structured()` |
| `psqn_generic()` | `psqn_optimize_generic()` |
| `psqn_bfgs()` | `psqn_bfgs_optimize()` |
| `psqn_aug_Lagrang()` | `psqn_aug_lagrang_structured()` |
| `psqn_aug_Lagrang_generic()` | `psqn_aug_lagrang_generic()` |
| `psqn_hess()` | `psqn_structured_hess()` |
| `psqn_generic_hess()` | `psqn_generic_hess()` |

The R Hessian routines return Eigen sparse matrices. The Fortran Hessian routines
currently return ordinary dense symmetric arrays; the numerical differentiation
and elementwise assembly are translated, while the Eigen sparse container is not.

## Compiler portability note

Version 0.1.1 avoids host association of an *optional* procedure dummy for the
constraint callback inside the optimizer's internal evaluation routine. This is
a source-level portability workaround for GNU Fortran versions that otherwise
report the constraint callback as having an implicit interface under
`-Werror=implicit-interface`.
