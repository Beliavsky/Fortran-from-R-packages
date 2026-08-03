# Testing

The permanent suite contains seven programs:

1. Exact Black-Scholes, default, dividend, delta, and vega benchmarks
2. Rate, volatility, survival, and price-linked hazard term structures
3. Dividend shifts and coupon present values
4. Grid construction, tridiagonal coefficients/solve, and European PDE pricing
5. American control-variate pricing and power-law default intensity
6. Zero-coupon, coupon, callable, and convertible bonds
7. Implied volatility, equivalent jump volatility, and finite-difference Greeks

Use:

```console
fpm test
```

The supplied scripts additionally compile all examples and the demo in checked
or optimized mode.
