# hermite-fortran

Modern Fortran/FPM translation of the computational code in the R package
`hermite` 1.2.1.

Implemented:

- exact generalized-Hermite PMF from its Poisson-component representation;
- upstream recurrence `int.hermite`;
- Edgeworth CDF approximation `edg`;
- Cornish-Fisher quantile approximation `cofi`;
- `dhermite`, `phermite`, `qhermite`, and `rhermite`;
- mean and variance helpers;
- raw-matrix generalized-Hermite regression with log/identity links;
- fixed or automatically selected integer order `m`;
- numerical Hessian/covariance, fitted values, LR test, p-value, and AIC.

The regression API uses a numeric design matrix rather than reproducing R
formula/model-frame/S3 machinery. No external runtime dependencies are
required.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

See `PORTING_NOTES.md` for exact-vs-approximation behavior and regression
details.
