# AEP-fortran

Modern Fortran/FPM translation of the computational core of the CRAN **AEP** package (version 0.1.4), for the asymmetric exponential-power distribution of Dongming and Zinde-Walsh (2009).

## API

- `daep(x, alpha, sigma, mu, epsilon [, log_pdf])`
- `paep(x, alpha, sigma, mu, epsilon [, log_p, lower_tail])`
- `qaep(u, alpha, sigma, mu, epsilon)`
- `raep(x, alpha, sigma, mu, epsilon)`
- `fitaep(x, fit [, starts, max_iter, tol])`
- `regaep(y, x, fit [, max_iter, tol])`

Scalar and rank-1 array overloads are provided for d/p/q. `raep` uses inverse-transform generation; this is distributionally equivalent to the R package's stable-mixture generator and is simpler and more robust in a standalone Fortran library.

`fitaep` and `regaep` port the package's iterative weighted estimation scheme, with native bounded scalar minimization, linear algebra, and observed-information calculations.

## Build

```sh
fpm test
fpm run --example demo_aep
```

No external dependencies are required.

## Numerical notes

The R package reports `n.p <- 3` in `fitaep` information criteria even though four AEP parameters are estimated. This port uses 4 parameters for AIC/CAIC/BIC/HQIC. The R regression code also computes its F statistic as a product by residual sum of squares; this port uses the standard regression F ratio. These corrections are documented rather than silently reproducing apparent source bugs.

Plotting, startup messages, R data-frame/formula presentation, and the bundled `plasma` dataset are not part of the numerical library.
