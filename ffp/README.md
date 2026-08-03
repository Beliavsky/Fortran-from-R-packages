# FFP Fortran

Modern Fortran/FPM translation of the computational core of the R package
**FFP: Fully Flexible Probabilities for Stress Testing and Portfolio Construction**.

## Included

- crisp and exponentially decaying flexible probabilities
- Gaussian-kernel probabilities
- effective number of scenarios and relative entropy
- weighted means, covariance, skewness, kurtosis, VaR and CVaR
- equality-constrained entropy pooling using the dual Newton method
- inequality/equality entropy pooling using exponentiated-gradient penalties
- least-information kernels and fitting probabilities to target moments
- double-decay covariance and probabilities
- weighted bootstrap scenario sampling
- constructors for mean, covariance, volatility, correlation and rank views
- half-life conversion

## Build

```console
fpm build
fpm test
fpm run
fpm run --example mean_view
```

The FPM version is numeric (`0.2.3`) for compatibility with FPM's manifest
parser.

See `PORTING.md`, `API.md`, and `TESTING.md` for details.
