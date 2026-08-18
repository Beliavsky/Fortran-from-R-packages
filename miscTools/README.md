# miscTools-fortran

Modern Fortran/FPM translation of the computational core of the R package
`miscTools` 0.6-30.

Implemented numerical functionality includes:

- derivatives of normal densities;
- column/row medians, including 3-D and 4-D array variants;
- matrix row/column insertion;
- symmetric-matrix construction and triangular/vector conversions;
- positive/negative semidefiniteness testing;
- quasiconcavity and quasiconvexity bordered-Hessian tests;
- R-squared;
- standard errors from covariance matrices;
- coefficient/t/p-value tables;
- string-name validation;
- simple observation/parameter counts.

The library is standalone and requires no runtime dependencies.

Plotting (`compPlot`, `histDens`), data-frame summary/digest reporting, R
attribute preservation, and generic fitted-model S3 infrastructure are
intentionally omitted.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The complete supplied R package is retained under `upstream/`.
