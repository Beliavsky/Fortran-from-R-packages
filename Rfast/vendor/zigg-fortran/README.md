# zigg-fortran

Modern Fortran 2018 translation of the computational core of the R package
`zigg` 0.0.2 by Dirk Eddelbuettel.

The package implements the Marsaglia-Tsang Ziggurat generator with the
Leong-Zhang et al. 32-bit correction used by upstream `zigg`.  It generates
standard normal, unit exponential, and uniform variates and preserves the
upstream 32-bit KISS/SHR3 state transitions.

## FPM

```sh
fpm build
fpm test
fpm run --example basic
```

## Module-level API

```fortran
use zigg

call zsetseed(12345)
x = zrnorm(1000)
y = zrexp(1000)
u = zrunif(1000)
```

The module-level routines use a saved default generator, mirroring the R
package's single static C++ generator.

## Independent generator objects

```fortran
type(ziggurat_rng) :: rng
real(dp) :: x(100)

call rng%set_seed(12345)
call rng%fill_normal(x)
```

Scalar methods `rnorm()`, `rexp()`, `runi()`, and `kiss()` are also available.
The four-word internal state can be saved and restored with `get_state` and
`set_state` (or `zgetstate`/`zsetstate` for the default generator).

## Compatibility

For seed 12345, the regression tests compare the first 20 draws of each
family directly with draws from the supplied upstream C++ header.  The
normal, exponential, uniform, and four-word internal state sequences match.

As upstream notes, this is a lightweight, fast historical generator with a
shorter period and potentially lower quality than modern general-purpose RNGs.
It should not be used for cryptography.

## License

GPL-2.0-or-later, following upstream `zigg`.  See `LICENSES.md` and `upstream/`.
