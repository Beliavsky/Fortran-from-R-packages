# ucminf-fortran

Modern free-form Fortran translation of the computational core of the R package
`ucminf` 1.2.3.

The package implements Hans Bruun Nielsen's UCMINF unconstrained nonlinear
optimizer: inverse-Hessian BFGS updating, soft line search, and trust-region-style
monitoring of the line-search input step. Analytic gradients and the original
forward/central finite-difference gradient choices are supported.

## Build with FPM

```sh
fpm build
fpm test
fpm run --example rosenbrock
```

## Minimal use

```fortran
use ucminf, only : dp, ucminf_result, ucminf_minimize

type(ucminf_result) :: result
real(dp) :: x0(2)

x0 = [2.0_dp, 0.5_dp]
call ucminf_minimize(x0, objective, result)
```

Supply an analytic gradient as the fourth argument when available. Otherwise the
optimizer uses the finite-difference method selected in `ucminf_options`.

For library integrations that need objective-specific data without globals or
compiler-generated nested-procedure trampolines, `ucminf_minimize_context` accepts
a polymorphic context object and callbacks of the form `fn(x, context)`. This is a
Fortran API extension; it does not change the translated UCMINF algorithm.

## Upstream provenance

The translation is based on `ucminf` 1.2.3. The original algorithm and Fortran
implementation are by Hans Bruun Nielsen, IMM/DTU. The R implementation is by
Stig Bousgaard Mortensen and subsequent maintainers/contributors are listed in
`upstream/DESCRIPTION` and the retained upstream sources.

The original fixed-form Fortran, R wrapper, C interface, DESCRIPTION, and README
are retained under `upstream/` for provenance and comparison.

## License

GPL-2.0-or-later, matching the upstream package's `GPL (>= 2)` license. See
`COPYING` and the upstream metadata.
