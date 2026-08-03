# Porting notes

## R objects

R lists, formula arguments, `...`, S3 `optimx`/`opm` objects, and data-frame
methods are represented by explicit Fortran derived types. The R replacement
method `coef<-` is not meaningful in a statically typed numerical library and
is the only exported namespace item not compiled.

## Solver mappings

The package coordinates many solvers from base R and optional packages. This
self-contained port implements the solver families directly:

- `Rvmmin`, BFGS, L-BFGS-B, `nvm`, and `nlminb` names use projected full-memory
  BFGS with an Armijo backtracking search.
- `Rcgmin`, CG, and `ncg` use projected Polak-Ribiere-plus nonlinear CG.
- `hjn` uses Hooke-Jeeves exploratory and pattern moves.
- `Nelder-Mead`, `anms`, and `nmsimplex` use one portable bounded simplex
  implementation.
- `snewton`, `snewtm`, `tn`, `tnbc`, Newton, SPG, and `ucminf` names use a
  safeguarded Hessian solve with diagonal regularization and line search.

Thus the high-level mathematical behavior is retained, but exact iteration
paths of base R, `nloptr`, `pracma`, or other optional backends are not claimed.

## Licensing

The upstream package declares GPL-2. This translation was written directly
from the optimx sources and standard published algorithms. Code from the
previous LGPL-3/GPL-3 translations was not copied into this GPL-2-only project.
