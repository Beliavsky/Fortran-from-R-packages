# BivGeo-fortran

Modern Fortran translation of the computational code in the R package `BivGeo` 2.1.1, which implements the Basu-Dhar bivariate geometric distribution.

The library is self-contained and uses only free-format Fortran 2018 source. No C, C++, R, or external numerical library is required.

## Main API

- `bivgeo_params`, `make_bivgeo_params`
- `dbivgeo1`, `dbivgeo2`
- `pbivgeo`, `sbivgeo`
- `cfbivgeo`, `covbivgeo`, `corbivgeo`
- `mean_bivgeo`, `variance_bivgeo`
- `mombivgeo`
- `rbivgeo1`, `rbivgeo2`
- `bivgeo_seed`

`rbivgeo1` uses exact discrete conditional inversion and `rbivgeo2` uses the Marshall-Olkin shock construction.

## Build

```text
fpm build
fpm test
fpm run --example example_bivgeo
```

## License

GPL-2.0-or-later, matching the upstream package's `GPL (>= 2)` declaration. See `LICENSE` and `UPSTREAM.md`.
