# RcppNumerical-fortran

Modern Fortran 2018 translation of the computational core of the R package
`RcppNumerical` 0.7-0.

The package provides callback-based numerical integration and optimization
without R, Rcpp, or Eigen:

- adaptive one-dimensional Gauss-Kronrod integration;
- all 12 public rules from 15 through 201 points;
- finite, semi-infinite, and doubly infinite intervals;
- deterministic multidimensional Cuhre cubature;
- unconstrained L-BFGS optimization;
- box-constrained L-BFGS-B optimization; and
- the package's `fastLR` logistic-regression example as a native Fortran API.

## Build with FPM

```text
fpm build
fpm test
fpm run --example integration_optimization
fpm run --example fast_lr_example
```

The optimizer implementations are bundled as local FPM dependencies under
`dependencies/lbfgs` and `dependencies/lbfgsb3`. Their package names match the
keys in the root `fpm.toml`.

## Basic API

```fortran
use rcppnumerical

type(integration_result_t) :: integral
type(optimization_result_t) :: optimum
real(dp) :: x(2)

call integrate_1d(my_scalar_function, 0.0_dp, 1.0_dp, integral)
call integrate_nd(my_vector_function, lower, upper, cubature)
call optim_lbfgs(my_objective_gradient, x, optimum)
call optim_lbfgsb(my_objective_gradient, x, lower, upper, optimum)
```

Callbacks may optionally receive polymorphic user data. See the examples and
`docs/API_MAP.md`.

## Scope

Rcpp classes, Eigen expression templates, R warnings, dynamic routine
registration, and the R-facing `fastLR()` list wrapper are not reproduced.
They are replaced by explicit Fortran procedure interfaces and derived result
types. No plotting code is present in the source package.

## Licensing

The translated package is distributed under GPL-2.0-or-later, matching the R
package. Portions retain their upstream MPL-2.0, LGPL-3.0, MIT, and BSD notices.
See `NOTICE.md`, `licenses/`, and the license files in each bundled dependency.
