# SkewHyperbolic-fortran

Modern Fortran/FPM implementation of the numerical core of the R package `SkewHyperbolic`.

## Implemented

- `dskewhyp`, vector density, and log density
- `ddskewhyp` density derivative
- `pskewhyp`, `qskewhyp`, and vector forms
- `rskewhyp` using the exact inverse-gamma/normal variance-mean mixture
- mean, variance, skewness, kurtosis, mode
- arbitrary integer moments about any centre
- numerical range calculation by density or tail probability
- maximum-likelihood fitting with positive `delta`/`nu` transforms
- automatic starting values

The parameter order is `(mu, delta, beta, nu)`, matching upstream.

## Build

```sh
fpm test
fpm run --example basic
```
