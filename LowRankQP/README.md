# LowRankQP-fortran

Modern Fortran/FPM translation of the computational code in R package
**LowRankQP 1.0.6** by John T. Ormerod and Matt P. Wand.

The package solves box-constrained equality-constrained convex quadratic
programs whose Hessian is either supplied directly as a square matrix or in
low-rank factor form.

## Algorithms translated

The C implementation was translated directly, including its primal-dual
predictor/corrector interior-point iteration and all four linear-system paths:

- full LU factorization;
- full Cholesky factorization;
- Sherman-Morrison-Woodbury (`SMW`);
- product-form Cholesky factorization (`PFCF`).

The original code called BLAS/LAPACK. This translation uses modern Fortran
array operations plus small self-contained LU and Cholesky solvers, so no
external numerical library is required.

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example quadprog_example
fpm run --example lowrank_example
```

## Minimal use

```fortran
use lowrankqp_kinds, only : dp
use lowrankqp, only : lowrankqp_options, lowrankqp_result, &
    solve_low_rank_qp, LRQP_PFCF

type(lowrankqp_options) :: opt
type(lowrankqp_result) :: res

opt%method = LRQP_PFCF
call solve_low_rank_qp(v, d, a, b, u, res, opt)
```

The Fortran interface expects `a(p,n)` with `a*alpha=b`. This is more natural
than the R/C boundary, where the R wrapper transposes its constraint matrix
before calling C.

## License

The upstream DESCRIPTION declares `GPL (>= 2)`. The translation is therefore
released under GPL-2.0-or-later. The complete supplied upstream package is
preserved under `original/LowRankQP-master/`.
