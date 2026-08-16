# polyaAeppli-fortran

Modern Fortran/FPM translation of `polyaAeppli` 2.0.2.

The Polya-Aeppli distribution is the geometric compound-Poisson distribution:
a Poisson number of independent shifted-geometric cluster sizes are summed.

Implemented:

- mass function and log mass;
- lower/upper CDF;
- log lower/upper tails;
- right-continuous integer quantiles;
- random generation;
- R-style vector parameter recycling;
- stable PMF, CDF, and upper-tail recurrence helpers;
- mean and variance helpers.

No external runtime dependencies are required.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The original R package is retained under `upstream/`.
