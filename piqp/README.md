# piqp-fortran

Modern Fortran/FPM translation of the computational interface of **piqp 0.6.2**
(PIQP-R), using the attached `Matrix-fortran` translation for CSC sparse-matrix
input.

The solver handles convex quadratic programs

```
minimize    0.5*x' P x + c' x
subject to  A x = b
            h_l <= G x <= h_u
            x_l <= x <= x_u
```

with positive-semidefinite `P`. Infinite lower/upper bounds are represented
with IEEE infinities.

## Build

```sh
fpm build
fpm test
```

## Dense example

```fortran
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
use piqp

real(dp) :: p(2,2), c(2), a(1,2), b(1), g(2,2)
real(dp) :: hu(2), xl(2), xu(2), inf
type(piqp_result_type) :: res

inf = ieee_value(0.0_dp, ieee_positive_inf)
p = reshape([6.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
c = [-1.0_dp,-4.0_dp]
a = reshape([1.0_dp,-2.0_dp],[1,2]); b = 0.0_dp
g = reshape([1.0_dp,-1.0_dp,0.0_dp,0.0_dp],[2,2])
hu = [1.0_dp,1.0_dp]
xl = [-inf,-1.0_dp]; xu = [inf,1.0_dp]

call solve_piqp(p,c,res,a,b,g,h_u=hu,x_l=xl,x_u=xu)
print *, res%x
```

The expected solution is approximately `(0.42857143, 0.21428571)`.

## Persistent model / updates

`piqp_model_type` mirrors the upstream setup/solve/update workflow:

```fortran
type(piqp_model_type) :: model
call model%setup(p,c,a,b,g,h_u=hu,x_l=xl,x_u=xu)
call model%solve()
call model%update(pmat=p_new, amat=a_new, h_u=hu_new, x_u=xu_new)
call model%solve()
```

The update API keeps dimensions fixed, as upstream PIQP does.

## Sparse input

The generic `solve_piqp` also accepts `matrix_sparse::csc_matrix` for `P`, `A`
and `G`. The first release converts CSC data to dense arrays and runs the same
proximal interior-point numerical core. This preserves sparse input semantics,
but not upstream PIQP's sparsity-preserving LDLT KKT factorization. See
`TRANSLATION_COVERAGE.md`.

## License

PIQP-derived `src/` files are BSD-2-Clause. The attached Matrix translation is
GPL-3.0-only, so the default combined FPM target is GPL-3.0-only. See
`NOTICE.md` and `LICENSES/`.
