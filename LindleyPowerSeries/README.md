# LindleyPowerSeries-fortran

Modern Fortran translation of the computational code in the R package
`LindleyPowerSeries` 1.0.1 by Saralees Nadarajah, Yuancheng Si, and Peihao
Wang.

The original package implements five members of the Lindley power-series
family described by Si and Nadarajah (2018). This translation preserves the
probability, density, hazard, quantile, and random-generation functionality
while replacing the R and `lamW` runtime dependencies with self-contained
Fortran 2018 code.

## Implemented distributions

- Lindley-geometric
- Lindley-logarithmic
- Lindley-negative-binomial
- Lindley-binomial
- Lindley-Poisson

For each family the module exports the same `p*`, `d*`, `h*`, `q*`, and `r*`
name used by the R package, with dots replaced by ordinary Fortran identifier
syntax where necessary. The scalar probability functions are `elemental`, so
they can be applied directly to conformable arrays.

The base Lindley CDF, PDF, and quantile are also public as `lindley_cdf`,
`lindley_pdf`, and `lindley_quantile`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example lps_demo
```

The project has no external dependencies.

## Example

```fortran
use lindley_power_series
implicit none

real(dp) :: x

x = qlindleypoisson(0.75_dp, 1.0_dp, 0.5_dp)
print *, x
print *, plindleypoisson(x, 1.0_dp, 0.5_dp)
```

## Numerical implementation

The inverse Lindley CDF uses a self-contained real `W_{-1}` Lambert-W
implementation with Halley iteration and a branch-point expansion. Stable
`log(1+x)` and `exp(x)-1` helpers are included because they are not standard
Fortran intrinsics. The Poisson and binomial power-series calculations use
log-space identities for large parameters, and the hazard functions use the
base Lindley survival/hazard directly to avoid `0/0` cancellation in far
tails.

See `TRANSLATION_NOTES.md` for the function map and two corrected upstream
formula defects.

## License

The upstream package declares `GPL (>= 2)`. This translation is distributed
under GPL-2.0-or-later and retains upstream authorship and source snapshots in
`upstream/` for provenance. See `LICENSE`.
