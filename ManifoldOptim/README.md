# ManifoldOptim-fortran 0.2.0

Modern Fortran/FPM translation of the computational interface of
**ManifoldOptim 1.0.2**, which embeds a ROPTLIB-derived Riemannian optimization
engine.

Version 0.2.0 closes the major computational architecture gaps in the initial
0.1.0 translation: constructed Stiefel retraction, sphere parameter-set
transport/locking behavior, LowRank quotient scaling, compact LRTRSR1,
ROPTLIB-style quasi-Newton secant scaling, multiple line-search algorithms and
a custom line-search hook.

The code is standalone Fortran 2018 and does not require R, Rcpp, Armadillo,
Eigen, BLAS, or LAPACK.

## Build

```text
fpm build
fpm test
fpm run --example sphere_example
fpm run --example brockett_stiefel
```

The package enables explicit-interface / implicit-typing safeguards in
`fpm.toml`.

## Basic use

```fortran
use manifoldoptim

type(manifold_domain) :: domain
type(solver_options) :: opt
type(solver_result) :: result

allocate(domain%component(1))
domain%component(1) = make_component(MANI_STIEFEL, n=10, p=3, param_set=2)

opt%tolerance = 1.0e-8_dp
call manifold_optimize(domain, x0, objective, gradient, 'LRBFGS', result, opt)
```

Callbacks use flat column-major arrays, matching the layout exposed by the R
package.

## Supported manifolds

- Euclidean
- Sphere
- Stiefel, including ParamSet 1 and constructed-retraction ParamSet 2
- Grassmann
- symmetric positive definite (SPD)
- LowRank (`U,D,V` factor layout with quotient-space tangent scaling)
- orthogonal group
- arbitrary product manifolds and multiplicities

## Supported solvers

`RSD`, `RCG`, `RNewton`, `RBFGS`, `LRBFGS`, `RBroydenFamily`, `RWRBFGS`,
`RTRSD`, `RTRNewton`, `RTRSR1`, and `LRTRSR1`.

If no analytical gradient is supplied, central numerical differences are used.
Newton methods without an analytical Hessian-vector callback use directional
finite differences of the Riemannian gradient.

See `API.md` and `TRANSLATION_COVERAGE.md` for details and the small remaining
implementation differences from ROPTLIB.

## License

The upstream package declares `GPL (>= 2)`. This translation is distributed
under GPL-2.0-or-later. See `LICENSE`, `licenses/GPL-2.0.txt`, and
`UPSTREAM_PROVENANCE.md`. The complete supplied source is retained in
`original/ManifoldOptim-master`.
