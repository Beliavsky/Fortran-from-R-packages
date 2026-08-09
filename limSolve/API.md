# API

All real arguments use `dp = kind(1.0d0)`.  Missing R matrices are represented
by zero-row arrays with the correct number of columns, e.g. `real(dp) :: E(0,n)`.

## Result types

`solve_result` contains `x`, `residual_norm`, `solution_norm`, `is_error`,
`status`, `numiter`, optional `unconstrained_solution`, optional `covariance`,
`rank_eq`, and `rank_app`.

`range_result` contains `range(:,1:2)` for min/max plus optional `central` and
`all_x`.

`sample_result` contains sampled `x`, optional null-space coordinates `q`,
probabilities `p`, jump lengths, acceptance ratio, and status.

## Main routines

```fortran
call solve_generalized(A, B, X [, tol] [, rank] [, status])
call nnls(A, B, result [, tol] [, max_iter])
call ldp(G, H, result [, tol] [, lower] [, upper])
call ldei(E, F, G, H, result [, tol] [, lower] [, upper])
call lsei(A, B, E, F, G, H, result [, tol] [, Wx] [, Wa] &
          [, lower] [, upper] [, fulloutput])
call linp(E, F, G, H, Cost, result [, ispos] [, int_vec] [, lower] [, upper])
```

Constraints follow the package convention `E*x = F`, `G*x >= H`.

For `lsei`, a scalar `Wx` requests unit-length column scaling and an `n`-vector
supplies explicit column scale factors, matching the Lawson-Hanson option keys
used by the R wrapper.  `Wa` scales rows of the approximate equations.

```fortran
r = resolution(A [, tol])
r = xranges(E,F,G,H [,ispos] [,central] [,full] [,lower] [,upper])
r = varranges(E,F,G,H,EqA [,EqB] [,ispos] [,lower] [,upper])
call varsample(X, EqA, Var [,EqB] [,status])
```

Linear systems:

```fortran
call solve_tridiag(dl, d, du, B, X, status)
call solve_banded(abd, nup, nlow, B, X, status [,full])
call solve_block(top, blocks, bot, B, overlap, X, status)
```

Sampling:

```fortran
call xsample(A,B,E,F,G,H,result, iter=..., outputlength=..., &
             type='mirror'|'rda'|'cda', jump=..., seed=...)
```
