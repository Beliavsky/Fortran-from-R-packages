# matrixNormal-fortran

Modern Fortran/FPM translation of the computational core of the R package `matrixNormal` 0.1.2.

## Included functionality

The public `matrixNormal` module provides:

* `identity_matrix`, `ones_matrix`
* `tr`, `vec`, `vech`
* `is_square_matrix`, `is_symmetric_matrix`
* `is_positive_definite`, `is_positive_semidefinite`
* `dmatnorm`
* `pmatnorm`
* `rmatnorm`
* `check_matnorm`
* the attached mvtnorm probability-control/result types and constructors (`genz_bretz`, `tvpack`, `miwa`)

The package uses the matrix-normal convention

```
vec(A) ~ N(vec(M), V kron U)
```

where `vec` stacks columns, exactly as R and Fortran store a matrix.

## Example

```fortran
program demo
  use matrixNormal
  implicit none
  real(dp) :: m(2,2), u(2,2), v(2,2), a(2,2)

  m = 0.0_dp
  u = reshape([1.0_dp,0.25_dp,0.25_dp,1.5_dp],[2,2])
  v = reshape([2.0_dp,0.30_dp,0.30_dp,1.0_dp],[2,2])

  a = rmatnorm(m,u,v,1234)
  print *, dmatnorm(a,m,u,v)
end program demo
```

For rectangle probabilities, `pmatnorm` returns the attached mvtnorm port's `probability_result`, so use `result%value`, `result%error`, and `result%inform`.

## Dependency

The user-supplied `mvtnorm-fortran` project is retained as a local FPM dependency under `vendor/mvtnorm-fortran`. It provides the multivariate Normal probability and simulation backend as well as numerical SPD/eigenvalue helpers.

## Deliberate source corrections

Two clear inconsistencies in the current R source are documented and corrected by default:

1. `pmatnorm` now uses the documented covariance `V kron U`; `legacy_covariance_order=.true.` reproduces the current R source's `U kron V` call.
2. `vech` performs actual lower-triangle half-vectorization instead of using lower-triangle numerical values as indices.

See `PORTING_NOTES.md` for details.

## Build

With FPM installed:

```text
fpm build
fpm test
fpm run --example basic
```

The project uses a relative local dependency path, so keep the included `vendor/mvtnorm-fortran` directory in place.

## Licensing

Top-level matrixNormal-derived Fortran sources are GPL-3.0-only, matching the upstream package's GPL-3 declaration. The bundled user-supplied mvtnorm dependency remains separately licensed GPL-2.0-only and retains its own complete license/notices. See `NOTICE.md` and `PORTING_NOTES.md`.
