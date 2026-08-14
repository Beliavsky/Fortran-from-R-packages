# mixsqp-fortran

Modern free-form Fortran/FPM translation of the computational code in the R
package `mixsqp` 0.3-54.

`mixsqp` solves the convex mixture-proportion problem

    minimize  -sum_i w_i log(sum_j L(i,j) x(j))
    subject to x(j) >= 0 and sum_j x(j) = 1.

The implementation follows the upstream mix-SQP algorithm: preliminary EM
updates, gradient/Hessian construction, an active-set quadratic subproblem,
and backtracking line search.

## Included computational functionality

- `fit_mixsqp`: high-level mixture-proportion solver.
- `mixobjective`: exported objective evaluator with R-compatible normalization
  of `x` and `w`.
- `mixem_update`: one EM update.
- `active_set_qp`: active-set quadratic solver used by mix-SQP.
- row likelihood and log-likelihood normalization.
- all-zero likelihood-column removal and solution reinsertion.
- numerical `eps` safeguard used by upstream mixsqp.
- optional low-rank SVD acceleration.
- final objective, gradient and Hessian.
- per-iteration objective, dual residual, nonzero count, step size, maximum
  iterate change, QP iteration count and line-search iteration count.
- `simulate_mix_data`, corresponding to upstream `simulatemixdata`.

All project source is `.f90` free source form. Implicit typing and implicit
external procedures are disabled.

## Example

```fortran
program example
  use mixsqp
  implicit none
  real(dp), allocatable :: y(:), s(:), L(:,:)
  type(mixsqp_result) :: fit
  type(mixsqp_control) :: ctl

  call set_seed(1)
  call simulate_mix_data(1000, 10, y, s, L)

  ctl = mixsqp_default_control()
  ctl%verbose = .false.
  call fit_mixsqp(L, fit, control=ctl)

  print *, trim(fit%status_message)
  print *, fit%x
end program
```

Build with FPM:

    fpm build
    fpm test
    fpm run --example basic

BLAS and LAPACK are required.

## Low-rank SVD note

The R package uses `irlba` to avoid computing a full SVD when a low-rank
approximation is requested. This standalone MIT-licensed translation uses
LAPACK `DGESDD` to obtain the same tolerance-based factorization semantics,
then retains only singular values greater than `tol_svd` (at least two, as in
upstream `tsvd`). This avoids bundling the separately licensed `irlba` port.
For very large matrices the separately translated `irlba-fortran` package can
be integrated later as an iterative backend without changing the SQP API.

## License

The upstream package is MIT licensed. See `LICENSE`, `LICENSE.upstream`, and
`UPSTREAM.md`.
