# bridgedist-fortran

Modern Fortran translation of the computational code in the R package
`bridgedist` 0.1.3 by Bruce Swihart.

The bridge distribution of Wang and Louis (2003) is implemented with the same
parameter `phi`, where `0 < phi < 1`.

## Implemented API

- `dbridge(x, phi, log_value)`
- `pbridge(q, phi, lower_tail, log_p)`
- `qbridge(p, phi, lower_tail, log_p)`
- `rbridge(x, phi)` for scalar `phi`
- `rbridge(x, phi(:))` with cyclic parameter recycling
- `dbridge_recycle`, `pbridge_recycle`, `qbridge_recycle`
- `bridge_mean`, `bridge_variance`

The scalar distribution functions are elemental, so normal Fortran scalar
expansion works naturally for conformable arrays. Explicit `*_recycle`
routines are supplied for users who need R-style cycling of unequal-length
argument vectors.

## Build

```sh
fpm build
fpm test
fpm run --example example_bridge
```

The project has no external dependencies.

## Numerical implementation

The formulas are algebraically the same as the R package, but the density and
CDF are evaluated in tail-stable forms. In particular, the CDF avoids direct
`exp(phi*q)` overflow, upper tails use distributional symmetry instead of
`1-F(q)`, and the log-density avoids `cosh(phi*x)` overflow.

## License

GPL-2.0-or-later, matching upstream `bridgedist`'s `GPL (>= 2)` declaration.
See `LICENSE`, `LICENSES.md`, and the retained source under `upstream/`.
