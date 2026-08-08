# trustOptim-fortran

Modern Fortran/FPM translation of the computational code in the R package
`trustOptim` 0.8.7.4 by Michael Braun.

The upstream package implements trust-region minimization using Steihaug
preconditioned conjugate gradients.  Three public optimization paths are
translated (with a generic `trust_optim` front end as well):

- `trust_optim_sr1` -- dense symmetric-rank-one quasi-Hessian updates;
- `trust_optim_bfgs` -- dense BFGS quasi-Hessian updates;
- `trust_optim_sparse` -- user-supplied sparse symmetric Hessians.

The sparse path stores only one triangle of the Hessian and uses sparse
Hessian-vector products in the trust-region CG iteration.  The modified
Cholesky preconditioner is formed as a dense matrix internally; this keeps the
package self-contained and does not densify the Hessian used by CG.

Also translated are the exported binary-choice demonstration routines
`binary.f`, `binary.grad`, and `binary.hess`, exposed as `binary_value`,
`binary_gradient`, and `binary_hessian`.

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example rosenbrock
fpm run --example sparse_quadratic
```

A strict GNU Fortran build script is available in `scripts/`.

## Compiler target

The source is standard Fortran 2018 and has been checked using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

No BLAS, LAPACK, Eigen, R, Rcpp, or RcppEigen dependency is required.

## License

The upstream package is licensed under MPL 2.0 or later; its computational
C++ headers specifically state MPL 2.0.  The translated source is distributed
under MPL-2.0 and retains upstream provenance.  The complete supplied upstream
archive is included under `original/trustOptim-master/`.
