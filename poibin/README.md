# poibin-fortran

Modern Fortran/FPM translation of the computational code in the R package
`poibin` 1.6 by Yili Hong.

The package implements the Poisson-binomial distribution with optional integer
multiplicity weights.  It provides PMF, CDF, quantile and random-generation
interfaces corresponding to `dpoibin`, `ppoibin`, `qpoibin` and `rpoibin`.

CDF methods retained from the R package:

- `DFT-CF`: exact inversion of the characteristic function
- `RF`: exact recursive formula
- `RNA`: refined normal approximation
- `NA`: normal approximation
- `PA`: Poisson approximation

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The library is self-contained and has no external dependencies.

## Minimal example

```fortran
use poibin, only : dp, dpoibin, ppoibin, qpoibin
real(dp) :: pp(3)
pp = [0.2_dp, 0.5_dp, 0.8_dp]
print *, dpoibin(2, pp)
print *, ppoibin(2, pp, 'DFT-CF')
print *, qpoibin(0.5_dp, pp)
```

## License

GPL-2.0-only, matching the upstream package. See `LICENSE` and
`TRANSLATION_NOTES.md` for provenance and implementation notes.
