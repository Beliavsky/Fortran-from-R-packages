# API

All real arguments use `real(dp)` where `dp = kind(1.0d0)`.

## Result types

`type(ls_result)` contains `x`, residuals, dual values, index ordering,
residual norm, objective, mode, iteration count, rank and `k`.

Status values are `LSEI_SUCCESS`, `LSEI_BAD_DIMENSIONS`,
`LSEI_ITERATION_LIMIT`, `LSEI_INFEASIBLE`, and `LSEI_NUMERICAL`.

## Solvers

- `nnls_solve(a,b,result[,max_iter])`
- `pnnls_solve(a,b,kfree,result[,sum_value,max_iter])`
- `ldp_solve(e,f,result[,tol])`
- `lsi_solve(a,b[,e,f],result[,lower,upper,tol])`
- `lsei_solve(a,b[,c,d,e,f],result[,lower,upper,tol])`
- `qp_solve(q,p,result[,c,d,e,f,lower,upper,tol])`
- `pnnqp_solve(q,p,kfree,result[,sum_value,tol])`
- `hfti_solve(a,b,result[,tol])`

Constraint conventions match R `lsei`: `C*x = d`, `E*x >= f`.

## Utilities

- `indx(x,v,ind)`
- `mat_maxs(x[,dim])`
