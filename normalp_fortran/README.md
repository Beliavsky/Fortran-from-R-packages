# normalp-fortran

Modern Fortran translation of the computational parts of the R package
`normalp` 0.7.2.1 by Angelo M. Mineo.

The package implements the exponential-power distribution (also called the
generalized error distribution) and the numerical estimation/regression
utilities supplied by `normalp`.

## Implemented API

| R routine | Fortran routine | Notes |
|---|---|---|
| `dnormp` | `dnormp` | scalar and rank-1 array interfaces |
| `pnormp` | `pnormp` | lower/upper tail and log probability |
| `qnormp` | `qnormp` | lower/upper tail and log probability |
| `rnormp` | `rnormp` | inverse-transform generator |
| `estimatep` | `estimatep` | inverse approximation and direct option |
| `paramp` | `paramp_fit` | fixed or estimated shape parameter |
| `kurtosis` | `kurtosis_p` | VI, B2 and Bp diagnostics |
| `lmp` | `lmp_fit` | Lp regression, fixed or estimated p |
| `simul.mp` | `simul_mp` | parameter-estimation simulation study |
| `simul.lmp` | `simul_lmp` | Lp-regression simulation study |

Plotting functions (`graphnp`, `qqnormp`, `qqlinep`, and S3 plot methods) and
R-specific print/summary/formula wrappers are intentionally not translated.
Their numerical content is available through the routines above.

## Numerical implementation

The R package delegates gamma probabilities and quantiles to R. This port is
self-contained and supplies a regularized incomplete-gamma implementation and
numerical gamma quantile inversion. No external numerical library is required.

The original `qnormp` source appears to apply `log()` rather than `exp()` when
`log.pr=TRUE`. This port uses the conventional R distribution-function meaning:
a supplied log probability is exponentiated before inversion.

## Build and test

```sh
fpm test
fpm run --example normalp_example
```

A strict GNU Fortran build was also tested with:

```sh
gfortran -std=f2018 -Wall -Wextra -fcheck=all
```

## License

The upstream DESCRIPTION declares `License: GPL`. This translation retains that
license declaration. See `NOTICE.md`, `UPSTREAM-DESCRIPTION`, and `COPYING`.
