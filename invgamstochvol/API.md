# API reference

## Kinds and status codes

```fortran
integer, parameter :: dp
integer, parameter :: invgam_success
integer, parameter :: invgam_invalid_argument
integer, parameter :: invgam_nonfinite_input
integer, parameter :: invgam_numerical_failure
```

## `ourgeo`

```fortran
value = ourgeo(a1, a2, b1, zstar [, niter] [, status])
```

Evaluates the first `niter` terms of

```text
2F1(a1, a2; b1; zstar).
```

The default is 500 terms. The translated package follows the upstream domain
used by the stochastic-volatility likelihood: `b1 > 0` and `0 <= zstar < 1`.

## `invgam_likelihood_result`

```fortran
type(invgam_likelihood_result)
   real(dp) :: total_loglik
   real(dp), allocatable :: loglik(:)
   real(dp), allocatable :: all_st(:)
   real(dp), allocatable :: all_ctil(:, :)
   real(dp), allocatable :: alogfac(:, :)
   real(dp), allocatable :: alogfac2(:)
   real(dp), allocatable :: alfac(:)
   integer :: nit
   integer :: niter
   integer :: status
   character(len=160) :: message
end type
```

The recursion arrays use zero lower bounds to preserve the indexing of the
original C++ implementation:

- `loglik(0:T-1)`
- `all_st(0:T)`
- `all_ctil(0:T-1,0:nit)`
- `alogfac(0:nit,0:max(nit,niter))`
- `alogfac2(0:max(nit,niter))`
- `alfac(0:max(nit,niter))`

## `lik_clo`

```fortran
call lik_clo(residuals, b2, nu, rho, result [, nit] [, niter] &
   [, nproc] [, nproc2])
```

Arguments:

- `residuals`: univariate residual series, length at least three
- `b2`: positive volatility-level parameter
- `nu`: positive inverse-gamma degrees-of-freedom parameter
- `rho`: persistence parameter satisfying `abs(rho) < 1`
- `nit`: outer series truncation, default 200
- `niter`: hypergeometric truncation, default 200
- `nproc`, `nproc2`: validated compatibility arguments; computation is serial

## `draw_k0`

Two typed forms are available.

```fortran
call draw_k0(result, nu, rho, b2, inverse_volatility &
   [, seed] [, status] [, nproc2])
```

or

```fortran
call draw_k0(all_st, all_ctil, alogfac, alogfac2, alfac, &
   nu, rho, b2, inverse_volatility [, seed] [, status] [, nproc2])
```

The result is a positive vector of length `T`, containing one exact posterior
draw of the inverse-volatility path. `seed` is an `integer(int64)` value.
