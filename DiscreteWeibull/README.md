# DiscreteWeibull-fortran

Modern Fortran/FPM translation of the computational code in the R package
`DiscreteWeibull` 1.1.

Implemented:

- type-I discrete Weibull on support `{1,2,...}` or `{0,1,...}`;
- type-III discrete Weibull;
- PMF, CDF, quantile, and random generation for both;
- type-III hazard;
- first/second moments, variance, and reciprocal moment where supplied;
- negative log-likelihood and moment-loss functions;
- `ML`, `M`, and `P` estimators for both families;
- type-I observed Fisher information and its inverse.

The supplied `Rsolnp-fortran` translation is vendored and used for the
type-I method-of-moments fit.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The library uses standard Fortran 2018. See `PORTING_NOTES.md` for numerical
details and the documented correction to the upstream type-III hazard.
