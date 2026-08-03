# Validation

The permanent tests cover:

1. matrix inversion, symmetric square roots, half-vectorization, polynomial
   products, pi weights, differencing, Corner tables, and ECCM tables;
2. VAR simulation and parameter recovery, forecasts, psi weights, impulse
   responses, FEVD, order selection, sparse/refined fitting, and Granger tests;
3. VARMA residuals, psi weights, covariance and conditional-likelihood fitting,
   plus VARX recovery, order selection, and multivariate regression;
4. cross-correlation and portmanteau/ARCH diagnostics, PCA/APCA, constrained
   factors, Stock-Watson forecasting, and Bayesian VAR;
5. EWMA fitting, DCC filtering/fitting, BEKK filtering/likelihood,
   modified-Cholesky and common-volatility procedures, and volatility
   diagnostics;
6. known-beta and Johansen-style VECM estimation and full/partial missing-value
   estimation.

Run checked and optimized validation with:

```text
./scripts/validate.sh
./scripts/validate_optimized.sh
```

The scripts compile and run all tests, the demonstration program, and every
example. Build products are written to temporary directories outside the source
tree.
