# API summary

## Nonlinear solver

### `nlsic_solve`

```fortran
call nlsic_solve(par0, m, residual, result, &
   jacobian, u, co, e, eco, control, mnorm)
```

Minimizes `||r(par)||_2` subject to

```text
u * par >= co
e * par  = eco
```

`jacobian`, `u/co`, `e/eco`, `control`, and `mnorm` are optional.

Important `nlsic_control` fields correspond directly to the R controls:

- `errx`
- `maxit`
- `btstart`, `btfrac`, `btdesc`, `btmaxit`, `btkmin`
- `rcond`
- `history`
- `adaptbt`
- `least_norm_step` (`sln` in R)
- `maxstep`
- `monotone`
- `reuse_jac`, `max_reuse`
- `report_ci`, `ci_p`

`nlsic_result` contains parameters, residuals, final Jacobian, last direction
and accepted step, covariance/CI information, iteration counters, convergence
status, and optional histories.

## Linear routines

All inequalities use the upstream convention `U*x >= co`.

```fortran
call ldp(u, co, result [, rcond])
call lsi(a, b, result [, u, co, rcond])
call lsi_ln(a, b, result [, u, co, rcond, mnorm, x0])
call lsie_ln(a, b, result [, u, co, e, ce, rcond, mnorm, x0])
call ls_ln(a, b, result [, rcond, mnorm, x0])
call ls_ln_multi(a, bmat, xmat, rank, status [, rcond, mnorm, x0])
call ls_ln_svd(a, b, result [, rcond])
call lsi_reg(a, b, result [, u, co, rcond, mnorm, x0])
call tls(a, b, x, status)
call nulla(v_or_matrix, basis, rank [, rcond])
call pnull(a, b, xp, basis, rank [, rcond, status])
call uplo_to_uco(lower, upper, u, co)
```

Use IEEE infinities in `lower`/`upper` for unspecified bounds.
