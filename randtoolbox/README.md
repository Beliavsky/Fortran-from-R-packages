# randtoolbox-fortran

Modern Fortran 2018 / FPM translation of the computational code in
`randtoolbox` 2.0.5.

The port is self-contained and does not require R.  It includes:

- general linear congruential generators, including exact modulo-2^64 arithmetic;
- SFMT for Mersenne exponents 607, 1279, 2281, 4253, 11213, 19937,
  44497, 86243, 132049, and 216091, including the upstream rotating parameter sets;
- MT19937 with the 2002 scalar and array initialization schemes and 32/53-bit uniforms;
- Knuth's TAOCP lagged-Fibonacci generator;
- all 17 WELL variants supplied by the `rngWELL` dependency;
- Torus, Halton, and 1111-dimensional Sobol low-discrepancy sequences;
- integer/bit conversion helpers and the first 100,000 prime numbers;
- gap, frequency, serial, poker, order, and collision-test computations;
- Stirling-number and permutation utilities used by the tests.

R's global `.Random.seed` integration (`set.generator`, `put.description`,
`get.description`, and the `.Call`/`.C` plumbing) is not reproduced.  Instead,
RNG state is held explicitly in Fortran derived types.  `trueRNG.R`, which
obtains numbers from an external web service, is also outside the numerical
library.  There is no plotting code in the package's computational port.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

An umbrella module is provided:

```fortran
use randtoolbox
```

The lower-level derived types (`sfmt_rng`, `mt19937_rng`, `well_rng`,
`congru_rng`, and `knuth_rng`) are public when exact state control or streaming
is needed.  Convenience array-returning procedures are provided by
`randtoolbox_pseudo`.

## Licensing

The upstream package is BSD-3-Clause.  The bundled MT19937 and SFMT notices,
Knuth public-domain notice, Sobol notice, and the separate WELL/rngWELL notices
are retained under `LICENSES/`.  See `UPSTREAM.md` for provenance.
